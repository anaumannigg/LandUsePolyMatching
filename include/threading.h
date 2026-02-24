#ifndef THREADING_H
#define THREADING_H

#include "command_line_parser.h"
#include "binary_io.h"
#include "cgal_includes.h"
#include "polygon_wh.h"
#include "shapefile_io_operations.h"
#include "localization.h"
#include "graph_computations.h"
#include "linear_program.h"

#include <string>
#include <condition_variable>


//define a struct for a decomposition task
struct DecompositionTask {
    //define min and max index of polys in set polys1 that should be queried for finding the intersected polys in polys2
    int min_idx,max_idx;
};

//define a struct for a decomposition from txt task
struct DecompositionFromTXTTask {
    std::pair<int, int> interval;
    int batch_size;
    std::string data_name;
};

//define a struct for a solving task
struct ConnectedSetTask {
    std::string data_name;
    double lambda;
    CommandLineOptions options;
    int task_id;
};

class Threading {
public:
    static void solvingWorkerTHREAD(
                  std::queue<ConnectedSetTask>& tasks,
                  std::mutex& queue_mutex,
                  std::condition_variable& cv,
                  std::atomic<bool>& done_flag,
                  std::atomic<int>& processed_counter);

    static void decompositionWorkerTHREAD(std::queue<DecompositionTask>& tasks,
                  std::mutex& queue_mutex,
                  std::condition_variable& cv,
                  std::atomic<bool>& done_flag,
                  const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2,
                  const Localization& rtree2, std::vector<Edge>& thread_computed_edges, std::atomic<int>& processed_counter);

    static void decompositionFromTXTWorkerTHREAD(std::queue<DecompositionFromTXTTask>& tasks,
                  std::mutex& queue_mutex,
                  std::condition_variable& cv,
                  std::atomic<bool>& done_flag,
                  const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2,
                  const std::string& filename, std::atomic<int>& processed_counter);

    static void statusTHREAD(std::atomic<int>& processed_counter, int num_sets);

    static void solveConnectedSet(ConnectedSetTask task, std::atomic<int>& processed_counter);

    static void decomposeIntoConnectedComponents(DecompositionFromTXTTask task, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const std::vector<Polygon_wh>& merged_polys,
                            const Localization& rtree1, const Localization& rtree2,
                            std::atomic<int>& processed_counter);

    static std::vector<Edge> computeIntersectionsAsEdges(DecompositionTask task, const std::vector<Polygon_wh>& polys1, const std::vector<Polygon_wh>& polys2, const Localization& rtree2, std::atomic<int>& processed_counter);

    static void decomposeIntoConnectedComponentsUsingTXT(DecompositionFromTXTTask task,
                                                 const std::vector<Polygon_wh>& polys1,
                                                 const std::vector<Polygon_wh>& polys2,
                                                 const std::string& filename,
                                                 std::atomic<int>& processed_counter);

    static Solution combineThreadSolutions(const std::string& data_name, int num_batches, int num_polys1, int num_polys2, std::vector<int>& set_ids,
                                        std::vector<std::pair<int,int>>& set_sizes,
                                        std::vector<std::vector<double>>& execution_times);

    static void cleanup_temp_folder(const std::string& folder_path);

private:
    //helper function to write thread solution to a file
    static void write_thread_results_to_file(const std::string& filename,
                                  const std::vector<int>& set_ids,
                                  const std::vector<std::pair<int,int>>& set_sizes,
                                  const std::vector<Solution>& sols,
                                  const std::vector<std::vector<double>>& execution_times);

    static void read_thread_results_from_file(const std::string& filename,
                                   std::vector<int>& set_ids,
                                   std::vector<std::pair<int,int>>& set_sizes,
                                   std::vector<Solution>& sols,
                                   std::vector<std::vector<double>>& execution_times);
};



#endif //THREADING_H
