#include "../include/cc_decomposition.h"

#include <condition_variable>

#include  "../include/union_find.h"
#include "../include/threading.h"

namespace CCDecomposition {
    void decomposeAndWriteToTXT(const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::string& filename, const int& num_threads) {
        //create R-Trees of the second dataset for effectively finding intersections
        Localization rtree2(polys2);

        //create tasks via fixed batch size
        int decomposition_batch_size = 10000;
        std::queue<DecompositionTask> decomp_tasks;
        for (int i=0; i < polys1.size(); i+=decomposition_batch_size) {
            decomp_tasks.push((DecompositionTask){i, std::min(i+decomposition_batch_size,(int)polys1.size())});
        }

        //create global storage for computed edges per thread (i.e., <int,int> pairs of poly indices per dataset)
        std::vector<std::vector<Edge>> edges_per_thread = std::vector<std::vector<Edge>>(num_threads);

        //start status thread
        std::atomic<int> processed_counter(0);
        std::thread status_thread(Threading::statusTHREAD, std::ref(processed_counter),polys1.size());

        //init mutexes
        std::mutex queue_mutex;
        std::condition_variable cv;
        std::atomic<bool> done_flag(false);

        //launch worker threads to compute edges
        std::vector<std::thread> decomp_threads;
        for(int t=0; t < num_threads; t++) {
            decomp_threads.emplace_back(Threading::decompositionWorkerTHREAD, std::ref(decomp_tasks),
                std::ref(queue_mutex), std::ref(cv), std::ref(done_flag), std::ref(polys1), std::ref(polys2), std::ref(rtree2), std::ref(edges_per_thread[t]),std::ref(processed_counter));
        }
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            done_flag = true;
        }
        cv.notify_all();

        //join threads
        for(int t=0; t < num_threads; t++) {
            decomp_threads[t].join();
        }

        //ensure the status thread finished
        processed_counter = polys1.size();
        status_thread.join();

        size_t num_polys1 = polys1.size();
        size_t num_polys2 = polys2.size();
        size_t total_num_polys = num_polys1 + num_polys2;

        // Create DSU
        UnionFind uf(total_num_polys);

        // Union all edges
        for (const auto& edges : edges_per_thread) {
            for (const auto& edge : edges) {
                int u = edge.first;
                int v = num_polys1 + edge.second;  // shift second dataset
                uf.unite(u, v);
            }
        }

        // Map root -> compact component id
        std::unordered_map<int, int> root_to_component;
        int component_counter = 0;

        std::vector<int> components(total_num_polys);

        for (size_t v = 0; v < total_num_polys; ++v) {
            int root = uf.find(v);

            if (root_to_component.find(root) == root_to_component.end()) {
                root_to_component[root] = component_counter++;
            }

            components[v] = root_to_component[root];
        }

        int num_computed_components = component_counter;

        std::vector<std::vector<int>> comp1(num_computed_components);
        std::vector<std::vector<int>> comp2(num_computed_components);

        for (size_t v = 0; v < total_num_polys; ++v) {
            int cid = components[v];

            if (v < num_polys1)
                comp1[cid].push_back(v);
            else
                comp2[cid].push_back(v - num_polys1);
        }

        writeDecompositionToTXT(num_computed_components, comp1, comp2, filename);
    }

    void writeDecompositionToTXT(const int& num_computed_components, const std::vector<std::vector<int>>& comp1, const std::vector<std::vector<int>> comp2, const std::string& filename) {
        //write to csv for future reference
        std::ofstream out(filename);

        for (int c = 0; c < num_computed_components; ++c) {

            if (comp1[c].empty() && comp2[c].empty())
                continue;

            // Write polys1 indices
            for (size_t i = 0; i < comp1[c].size(); ++i) {
                out << comp1[c][i];
                if (i + 1 < comp1[c].size())
                    out << ",";
            }

            out << " ";

            // Write polys2 indices
            for (size_t i = 0; i < comp2[c].size(); ++i) {
                out << comp2[c][i];
                if (i + 1 < comp2[c].size())
                    out << ",";
            }

            out << "\n";
        }
    }
}
