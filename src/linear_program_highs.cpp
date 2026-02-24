#include "../include/linear_program.h"
#include "Highs.h"
#include "interfaces/highs_c_api.h"

// ---- HiGHS VERSIONS ----
// HiGHS uses a large number for infinity if not explicitly set
static constexpr double kInf = 1e100;

void LinearProgram::solveILP(const Graph& g,
                            int num_polys1,
                            int num_polys2,
                            Solution* solution) {
  try {
    Highs highs;
    highs.setOptionValue("output_flag", false);


    std::vector<HighsInt> var;     // column indices
    std::vector<double> weight;

    std::vector<std::vector<HighsInt>> var_adj1_v(g.vertex_set().size());
    std::vector<std::vector<HighsInt>> var_adj2_v(g.vertex_set().size());

    std::vector<int> sources, targets;

    Graph::vertex_iterator v, vend;
    for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
      if (!g[*v].referenced_map) {
        for (auto e = out_edges(*v, g); e.first != e.second; ++e.first) {
          const double w = g[*e.first].weight;

          // --- Add binary var (version-compatible) ---
          const HighsInt col = highs.getNumCol();  // index of the next column to be added

          HighsStatus st = highs.addVar(0.0, 1.0); // returns status in your version
          if (st != HighsStatus::kOk) {
            std::cout << "HiGHS addVar failed\n";
            return;
          }

          st = highs.changeColCost(col, w);
          if (st != HighsStatus::kOk) {
            std::cout << "HiGHS changeColCost failed\n";
            return;
          }

          st = highs.changeColIntegrality(col, HighsVarType::kInteger);
          if (st != HighsStatus::kOk) {
            std::cout << "HiGHS changeColIntegrality failed\n";
            return;
          }
          // -----------------------------------------

          var.push_back(col);
          weight.push_back(w);

          var_adj1_v[*v].push_back(col);
          var_adj2_v[e.first->m_target].push_back(col);

          sources.push_back(*v);
          targets.push_back(e.first->m_target);
        }
      }
    }

    highs.changeObjectiveSense(ObjSense::kMaximize);

    // Constraints: each polygon picked at most once
    for (bool map_switch : {false, true}) {
      const int num_polys = !map_switch ? num_polys1 : num_polys2;

      for (int i = 0; i < num_polys; i++) {
        std::vector<int> referring_vertices;
        for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
          if (g[*v].referenced_map == map_switch &&
              std::find(g[*v].referenced_polys.begin(),
                        g[*v].referenced_polys.end(),
                        i) != g[*v].referenced_polys.end()) {
            referring_vertices.push_back(*v);
          }
        }

        if (referring_vertices.empty()) continue;

        const auto& adjacent_edges = !map_switch ? var_adj1_v : var_adj2_v;

        std::vector<HighsInt> idx;
        std::vector<double> val;
        for (int rv : referring_vertices) {
          for (HighsInt adj_var : adjacent_edges[rv]) {
            idx.push_back(adj_var);
            val.push_back(1.0);
          }
        }
        if (idx.empty()) continue;

        HighsStatus st = highs.addRow(-kInf, 1.0, (HighsInt)idx.size(), idx.data(), val.data());
        if (st != HighsStatus::kOk) {
          std::cout << "HiGHS addRow failed\n";
          return;
        }
      }
    }

    // Solve
    HighsStatus run_status = highs.run();
    if (run_status != HighsStatus::kOk) {
      std::cout << "HiGHS run() failed\n";
      return;
    }

    HighsModelStatus ms = highs.getModelStatus();
    if (!(ms == HighsModelStatus::kOptimal)) {
      std::cout << "HiGHS did not find a feasible solution. ModelStatus=" << (int)ms << "\n";
      return;
    }

    const HighsSolution& sol = highs.getSolution();
    if (sol.col_value.empty()) {
      std::cout << "HiGHS returned empty solution\n";
      return;
    }

    // Retrieve solution
    for (size_t i = 0; i < var.size(); i++) {
      const double x = sol.col_value[var[i]];
      if (x > 0.5) {
        double w = 0.0;
        bool weight_found = false;

        for (auto e = out_edges(sources[i], g); e.first != e.second; ++e.first) {
          if ((int)e.first->m_target == targets[i]) {
            w = g[*e.first].weight;
            weight_found = true;
            break;
          }
        }
        assert(weight_found && "did not find weight when iterating through edges!");

        solution->addMatch(g[sources[i]].referenced_polys,
                           g[targets[i]].referenced_polys,
                           w);
      }
    }

  } catch (const std::exception& e) {
    std::cout << "Exception during optimization (HiGHS): " << e.what() << "\n";
  } catch (...) {
    std::cout << "Unknown exception during optimization (HiGHS)\n";
  }
}

void LinearProgram::solveILP_trees(TreeConstrainedCandidateGraph& cg_tree,
                                  int num_polys1,
                                  int num_polys2,
                                  Solution* solution) {
  Graph& g = cg_tree.get_graph();

  // make sure the graph is directed upwards from the leafs to ensure an efficient collection of constraints
  cg_tree.setTreeDirection(true);

  try {
    Highs highs;
    highs.setOptionValue("output_flag", false);


    // variable column indices (one per edge variable)
    std::vector<HighsInt> var;

    // adjacency lists of variables per vertex (split by map side)
    std::vector<std::vector<HighsInt>> var_adj1_v(g.vertex_set().size());
    std::vector<std::vector<HighsInt>> var_adj2_v(g.vertex_set().size());

    // remember edge endpoints for solution retrieval
    std::vector<int> sources, targets;

    // ----- Create binary variables per outgoing edge of left side -----
    Graph::vertex_iterator v, vend;
    for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
      if (!g[*v].referenced_map) {
        for (auto e = out_edges(*v, g); e.first != e.second; ++e.first) {
          const double w = g[*e.first].weight;

          const HighsInt col = highs.getNumCol();

          HighsStatus st = highs.addVar(0.0, 1.0);  // bounds for binary
          if (st != HighsStatus::kOk) {
            std::cout << "HiGHS addVar failed\n";
            return;
          }

          st = highs.changeColCost(col, w);
          if (st != HighsStatus::kOk) {
            std::cout << "HiGHS changeColCost failed\n";
            return;
          }

          st = highs.changeColIntegrality(col, HighsVarType::kInteger);
          if (st != HighsStatus::kOk) {
            std::cout << "HiGHS changeColIntegrality failed\n";
            return;
          }

          var.push_back(col);

          var_adj1_v[*v].push_back(col);
          var_adj2_v[e.first->m_target].push_back(col);

          sources.push_back(*v);
          targets.push_back(e.first->m_target);
        }
      }
    }

    // Maximize objective
    highs.changeObjectiveSense(ObjSense::kMaximize);

    // ----- Add legality constraints using tree paths -----
    for (bool map_switch : {false, true}) {
      const int num_polys = !map_switch ? num_polys1 : num_polys2;

      for (int i = 0; i < num_polys; i++) {
        // collect all vertices referring to polygon i via path-to-root
        std::vector<Vertex> referring_vertices =
            cg_tree.getPathToRoot(map_switch, i, PathtoRootMode::START_AT_POLYGON_ID);

        if (referring_vertices.empty()) continue;

        const auto& adjacent_edges = !map_switch ? var_adj1_v : var_adj2_v;

        std::vector<HighsInt> idx;
        std::vector<double> val;

        for (const auto& rv : referring_vertices) {
          for (HighsInt adj_var : adjacent_edges[rv]) {
            idx.push_back(adj_var);
            val.push_back(1.0);
          }
        }

        if (idx.empty()) continue;

        // sum <= 1
        HighsStatus st = highs.addRow(-kInf, 1.0, (HighsInt)idx.size(), idx.data(), val.data());
        if (st != HighsStatus::kOk) {
          std::cout << "HiGHS addRow failed\n";
          return;
        }
      }
    }

    // ----- Solve -----
    HighsStatus run_status = highs.run();
    if (run_status != HighsStatus::kOk) {
      std::cout << "HiGHS run() failed\n";
      return;
    }

    HighsModelStatus ms = highs.getModelStatus();
    if (!(ms == HighsModelStatus::kOptimal)) {
      std::cout << "HiGHS did not find a feasible solution. ModelStatus=" << (int)ms << "\n";
      return;
    }

    const HighsSolution& sol = highs.getSolution();
    if (sol.col_value.empty()) {
      std::cout << "HiGHS returned empty solution\n";
      return;
    }

    // ----- Retrieve solution -----
    for (size_t i = 0; i < var.size(); i++) {
      const double x = sol.col_value[var[i]];

      // Use tolerance (don't compare doubles to exactly 1.0)
      if (x > 0.5) {
        double w = 0.0;
        bool weight_found = false;

        for (auto e = out_edges(sources[i], g); e.first != e.second; ++e.first) {
          if ((int)e.first->m_target == targets[i]) {
            w = g[*e.first].weight;
            weight_found = true;
            break;
          }
        }

        assert(weight_found && "did not find weight when iterating through edges!");

        solution->addMatch(g[sources[i]].referenced_polys,
                           g[targets[i]].referenced_polys,
                           w);
      }
    }

  } catch (const std::exception& e) {
    std::cout << "Exception during optimization (HiGHS): " << e.what() << "\n";
  } catch (...) {
    std::cout << "Unknown exception during optimization (HiGHS)\n";
  }
}



std::vector<double> solveLP_trees_fractionally(TreeConstrainedCandidateGraph& cg_tree,
                                               int num_polys1,
                                               int num_polys2) {
  Graph& g = cg_tree.get_graph();

  // make sure the graph is directed upwards from the leafs to ensure an efficient collection of constraints
  cg_tree.setTreeDirection(true);

  std::vector<double> fractional_values(cg_tree.getMaxEdgeID(), -1.0);

  try {
    Highs highs;
    highs.setOptionValue("output_flag", false);


    // variable (column) per edge
    std::vector<HighsInt> var;
    std::vector<double> weight;

    // adjacency lists: variables incident to a vertex (split by map side)
    std::vector<std::vector<HighsInt>> var_adj1_v(g.vertex_set().size());
    std::vector<std::vector<HighsInt>> var_adj2_v(g.vertex_set().size());

    // remember edge ids for solution retrieval
    std::vector<int> edge_ids;
    edge_ids.reserve(num_edges(g));

    // ----- Create variables per edge -----
    for (auto ei = edges(g); ei.first != ei.second; ++ei.first) {
      auto e = *ei.first;

      edge_ids.push_back(g[e].id);

      const double w = g[e].weight;

      // Version-compatible addVar: returns status, not the column index
      const HighsInt col = highs.getNumCol();

      HighsStatus st = highs.addVar(0.0, 1.0);  // continuous in [0,1]
      if (st != HighsStatus::kOk) {
        std::cout << "HiGHS addVar failed\n";
        return fractional_values;
      }

      st = highs.changeColCost(col, w);
      if (st != HighsStatus::kOk) {
        std::cout << "HiGHS changeColCost failed\n";
        return fractional_values;
      }

      // (No integrality call here: this is an LP, continuous vars)
      var.push_back(col);
      weight.push_back(w);

      // ensure source is in map1 and target is in map2
      Vertex source = !g[e.m_source].referenced_map ? e.m_source : e.m_target;
      Vertex target =  g[e.m_source].referenced_map ? e.m_source : e.m_target;

      var_adj1_v[source].push_back(col);
      var_adj2_v[target].push_back(col);
    }

    // maximize objective
    highs.changeObjectiveSense(ObjSense::kMaximize);

    // ----- Add legality constraints -----
    for (bool map_switch : {false, true}) {
      const int num_polys = !map_switch ? num_polys1 : num_polys2;

      for (int i = 0; i < num_polys; i++) {
        // vertices on the path to root for polygon i in the chosen map
        std::vector<Vertex> referring_vertices =
            cg_tree.getPathToRoot(map_switch, i, PathtoRootMode::START_AT_POLYGON_ID);

        if (referring_vertices.empty()) continue;

        const auto& adjacent_edges = !map_switch ? var_adj1_v : var_adj2_v;

        std::vector<HighsInt> idx;
        std::vector<double> val;

        for (const auto& rv : referring_vertices) {
          for (HighsInt adj_var : adjacent_edges[rv]) {
            idx.push_back(adj_var);
            val.push_back(1.0);
          }
        }

        if (idx.empty()) continue;

        // sum <= 1
        HighsStatus st = highs.addRow(-kInf, 1.0, (HighsInt)idx.size(), idx.data(), val.data());
        if (st != HighsStatus::kOk) {
          std::cout << "HiGHS addRow failed\n";
          return fractional_values;
        }
      }
    }

    // ----- Solve -----
    HighsStatus run_status = highs.run();
    if (run_status != HighsStatus::kOk) {
      std::cout << "HiGHS run() failed\n";
      return fractional_values;
    }

    HighsModelStatus ms = highs.getModelStatus();
    if (!(ms == HighsModelStatus::kOptimal)) {
      std::cout << "HiGHS did not find a feasible solution. ModelStatus=" << (int)ms << "\n";
      return fractional_values;
    }

    const HighsSolution& sol = highs.getSolution();
    if (sol.col_value.empty()) {
      std::cout << "HiGHS returned empty solution\n";
      return fractional_values;
    }

    // ----- Retrieve fractional solution -----
    for (size_t i = 0; i < var.size(); i++) {
      const HighsInt col = var[i];
      const double x = sol.col_value[col];

      const int id = edge_ids[i];
      if (0 <= id && id < (int)fractional_values.size()) {
        fractional_values[id] = x;
      }
    }

  } catch (const std::exception& e) {
    std::cout << "Exception during optimization (HiGHS): " << e.what() << "\n";
  } catch (...) {
    std::cout << "Unknown exception during optimization (HiGHS)\n";
  }

  return fractional_values;
}

//the matching algorithm proposed in 'On Tree Constrained Matchings and Generalizations' (Algorithm 1)
std::vector<int> CanzarMatching(TreeConstrainedCandidateGraph cg_tree, int num_polys1, int num_polys2) {
    //get optimal fractional solution
    auto X = solveLP_trees_fractionally(cg_tree, num_polys1, num_polys2);

    //first check if the solution is integer already, as then it can be returned right away
    bool is_integer=true;
    for(const auto& x : X) {
        if(x!= 0.0 && x!=1.0 && x != -1.0) {
            is_integer = false;
            break;
        }
    }
    if(is_integer) {
        std::vector<int> M;
        for(int e=0; e < X.size(); e++) {
            if(X[e] == 1.0) M.push_back(e);
        }
        return M;
    }

    //collect set of edges with fractional value 0
    std::vector<int> F_0;
    for(int x_id = 0; x_id < X.size(); x_id++) if(X[x_id] == 0.0) F_0.push_back(x_id);

    //check if set of 0-edges is not empty
    if(!F_0.empty()){
        //try to solve the matching problem again without the 0-weighted edges
        //create a deep copy of the tree
        TreeConstrainedCandidateGraph cg_tree_reduced(cg_tree);

        //collect sources and targets of all 0-edges
        std::vector<int> sources,targets;

        for (auto ei = edges(cg_tree_reduced.get_graph()); ei.first != ei.second; ++ei.first) {
            auto e = *ei.first;
            int e_id = cg_tree_reduced.get_graph()[e].id;

            if(std::binary_search(F_0.begin(),F_0.end(),e_id)) {
                //edge should be remembered for deletion
                sources.push_back(e.m_source);
                targets.push_back(e.m_target);
            }
        }

        //remove all 0-edges from the graph
        for(int e = 0; e<sources.size();e++) {
            cg_tree_reduced.delete_edge(sources[e],targets[e]);
        }

        //cout << "removed edges and recursed on graph with " << boost::num_edges(cg_tree_reduced.get_graph()) << " edges" << endl;
        //recurse
        auto M = CanzarMatching(cg_tree_reduced, num_polys1,num_polys2);
        return M;

    }

    //check if there exists an edge, such that all its conflicting edges are in sum < alpha = 3 in fractional value
    //conflicts are always the paths to the root
    int candidate_id = -1;
    double w_e = 0.0;
    std::vector<Graph::edge_descriptor> conflicting_edges;
    std::vector<int> conflicting_edge_ids;

    int dbg_edge_id = -1;
    for (auto ei = edges(cg_tree.get_graph()); ei.first != ei.second; ++ei.first) {
        dbg_edge_id++;
        auto e = *ei.first;

        //skip edges with decision variable value 1 as those have no conflict and thus no recursion is needed
        if(X[cg_tree.get_graph()[e].id] == 1.0) continue;

        //make sure source is 1 and target is 2
        Vertex source = !cg_tree.get_graph()[e.m_source].referenced_map ? e.m_source : e.m_target;
        Vertex target = cg_tree.get_graph()[e.m_source].referenced_map ? e.m_source : e.m_target;

        //collect paths to root
        std::vector<Vertex> root_path1 = cg_tree.getPathToRoot(0,source,PathtoRootMode::START_AT_VERTEX);
        std::vector<Vertex> subtree1 = cg_tree.getSubtree(0,source);
        std::vector<Vertex> root_path2 = cg_tree.getPathToRoot(1,target,PathtoRootMode::START_AT_VERTEX);
        std::vector<Vertex> subtree2 = cg_tree.getSubtree(1,target);

        //also collect all adjacent vertices of source and target as those pose conflicts as well
        auto adj_source = boost::adjacent_vertices(source,cg_tree.get_graph());
        std::vector<Vertex> adj2; for (auto ai = adj_source.first; ai != adj_source.second; ai++) {adj2.push_back(*ai);}
        auto adj_target = boost::adjacent_vertices(target,cg_tree.get_graph());
        std::vector<Vertex> adj1; for (auto ai = adj_target.first; ai != adj_target.second; ai++) {adj1.push_back(*ai);}



        std::vector<Vertex> vertices1 = root_path1;
        vertices1.insert(vertices1.end(),subtree1.begin(),subtree1.end());
        vertices1.insert(vertices1.end(),adj1.begin(),adj1.end());
        std::vector<Vertex> vertices2 = root_path2;
        vertices2.insert(vertices2.end(),subtree2.begin(),subtree2.end());
        vertices2.insert(vertices2.end(),adj2.begin(),adj2.end());

        //sort and make unique
        std::sort(vertices1.begin(),vertices1.end());
        vertices1.erase(std::unique(vertices1.begin(),vertices1.end()), vertices1.end());
        std::sort(vertices2.begin(),vertices2.end());
        vertices2.erase(std::unique(vertices2.begin(),vertices2.end()), vertices2.end());


        double conflict_sum = 0.0;

        //get all conflicts, meaning the edges adjacent to an 1 and/or 2 vertices out of the compiled set
        for (const auto& o : vertices1) {
            auto adj_edges = boost::out_edges(o,cg_tree.get_graph());
            for (auto ei = adj_edges.first; ei != adj_edges.second; ei++) {
                auto e_n = *ei;
                //skip the edge itself
                if (e_n == e) continue;
                conflict_sum += cg_tree.get_graph()[e_n].weight;
                conflicting_edges.push_back(e_n);
                conflicting_edge_ids.push_back(cg_tree.get_graph()[e_n].id);

            }
        }

        //note that we now only want to consider the edges that are not already considered (which may happen for edges incident to one of the vertices1)
        for (const auto& a : vertices2) {
            auto adj_edges = boost::out_edges(a,cg_tree.get_graph());
            for (auto ei = adj_edges.first; ei != adj_edges.second; ei++) {
                auto e_n = *ei;
                //skip the edge itself
                if (e_n == e) continue;
                //only add conflict if not considered yet
                if (std::find(conflicting_edges.begin(),conflicting_edges.end(),e_n) == conflicting_edges.end()) {
                    conflict_sum += cg_tree.get_graph()[e_n].weight;
                    conflicting_edges.push_back(e_n);
                    conflicting_edge_ids.push_back(cg_tree.get_graph()[e_n].id);
                }

            }
        }

        //if conflict sum is less than alpha = 3, a candidate is found
        if(conflict_sum < 3.0) {
            //cout << "setting candidate: " << e << endl;
            candidate_id = cg_tree.get_graph()[e].id;
            w_e = cg_tree.get_graph()[e].weight;
            break;
        }


    }

    //check if candidate has been found
    if(candidate_id > -1) {

        //modify weights of all edges in conflict with e
        for(const auto& n_e : conflicting_edges) {
            cg_tree.get_graph()[n_e].weight -= w_e;
        }

        //recurse
        auto M = CanzarMatching(cg_tree,num_polys1,num_polys2);

        //cout << "after recursion: " << boost::num_edges(cg_tree.get_graph()) << " edges in graph" << endl;

        //check if none of the edges in the conflicts are selected anymore
        std::sort(conflicting_edge_ids.begin(),conflicting_edge_ids.end());
        std::vector<int> edges_still_selected;
        std::set_intersection(conflicting_edge_ids.begin(),conflicting_edge_ids.end(),M.begin(),M.end(),std::back_inserter(edges_still_selected));
        if(edges_still_selected.empty()) {
            M.push_back(candidate_id);
            std::sort(M.begin(),M.end());
            return M;
        }

    }
    //arrived in base case
    //here, an empty matching should be returned
    return {};
}

void LinearProgram::solveViaCanzar(TreeConstrainedCandidateGraph& cg_tree, int num_polys1, int num_polys2, Solution& solution) {
    Graph& g = cg_tree.get_graph();
    //make sure the edges in the tree point upwards to the root
    cg_tree.setTreeDirection(true);

    //approximate Matching using Canzar's Algorithm
    auto M = CanzarMatching( cg_tree,num_polys1,num_polys2);

    //collect all found matches and add them to the solution
    // if no edges were selected, return
    if(M.size() == 0) return;
    int m_id = 0;
    //iterate over all edges in G
    for (auto ei = edges(cg_tree.get_graph()); ei.first != ei.second; ++ei.first) {
        auto e = *ei.first;
        int e_id = g[e].id;
        //M is sorted
        while(m_id < M.size()-1 && M[m_id] < e_id) m_id++;

        //check if edge is in matching
        if(e_id == M[m_id]) {
            //collect sets of polys to add to the solution
            bool source_map = g[e.m_source].referenced_map;
            std::vector<int> polys1 = !source_map? g[e.m_source].referenced_polys : g[e.m_target].referenced_polys;
            std::vector<int> polys2 = source_map? g[e.m_source].referenced_polys : g[e.m_target].referenced_polys;

            //add to solution
            solution.addMatch(polys1,polys2,g[e].weight);
        }
    }
}