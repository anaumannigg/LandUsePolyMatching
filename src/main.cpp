#include "../include/threading.h"

#include <exception>
#include <filesystem>
namespace fs = std::filesystem;

//measuring the peak memory used by the process
size_t get_peak_memory_kb() {
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.substr(0, 6) == "VmHWM:") { // VmHWM = High Water Mark (peak resident set size)
            std::istringstream iss(line);
            std::string key, value, unit;
            iss >> key >> value >> unit;
            return std::stoul(value); // in kilobytes
        }
    }
    return 0;
}

//takes two maps and a third map that contains the merged version of both maps with a union of all polygons with a common boundary
//decomposes the maps w.r.t. the merged map into separated subsets and saves them in a subfolder 'tmp' of the current folder
void run_NtoM_with_PreDecomposition(const CommandLineOptions& options) {
    //extract options
    std::string data_name = options.dataset_name;
    int num_threads = options.num_threads;
    SOLUTION_MODE sol_mode = options.mode;
    TREE_BUILD_MODE tree_mode = options.tree_mode;

    //retreive paths
    std::string data_descr1 = options.data_descr1;
    std::string data_descr2 = options.data_descr2;
    std::string file1 = "../input/" + data_name + "/" + data_name + "_" + data_descr1;
    std::string file2 = "../input/" + data_name + "/" + data_name + "_" + data_descr2;

    std::string folderTMP = "../input/" + data_name + "/tmp/";
    //create folder of it does not exist yet
    if (!fs::exists(folderTMP)) fs::create_directory(folderTMP);

    //overall timekeeping
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //read polygons from shapefile (if shapefile does not exist, try to read from gpkg
    std::vector <Polygon_wh> polys1, polys2, merged_polys;
    int num_invalid1 = 0, num_invalid2 = 0;
    try {
        num_invalid1 = ReadGeoPackage(file1 + ".gpkg", polys1);
        num_invalid2 = ReadGeoPackage(file2 + ".gpkg", polys2);
    } catch (std::exception& e) {
        std::cerr << "Error reading input files: " << e.what() << endl;
        return;
    }

    cout << "read " << polys1.size() << " polygons from file 1. (" << num_invalid1 << " ignored due to invalid geometries)" << endl;
    cout << "read " << polys2.size() << " polygons from file 2. (" << num_invalid2 << " ignored due to invalid geometries)" << endl;


    //if a file with the decomposition is already provided, this is the simplest way of getting the decomposition, prefer it
    std::string fileDECOMPOSITION = "../input/" + data_name + "/connected_components.txt";
    //if a file with the merged information is already provided, load it, else compute the union and save it for
    //later use
    std::string fileMERGED = "../input/" + data_name + "/" + data_name + "_merged";

    //store info about decomposed input
    double timing_decomp = 0.0;
    int num_batches = 0;
    int num_components = 0;
    if (fs::exists(fileDECOMPOSITION)) {
        std::chrono::steady_clock::time_point decomp_begin = std::chrono::steady_clock::now();

        //first count the lines (being the components)
        std::string unused;
        std::ifstream input_file(fileDECOMPOSITION);
        while (std::getline(input_file,unused)) ++num_components;

        //set sets per task
        int sets_per_task = options.batch_size;
        num_batches = num_components / sets_per_task;
        //there might be an additional smaller batch to consider all polygons
        bool last_batch_is_smaller = false;
        if (num_components % sets_per_task != 0) {
            num_batches++;
            last_batch_is_smaller = true;
        }

        std::vector<std::pair<int,int>> task_set_intervals;
        for(int t=0; t < num_batches;t++) {
            if (last_batch_is_smaller && t == num_batches - 1) {
                task_set_intervals.emplace_back((num_batches-1) * sets_per_task, num_components - 1);
                break;
            }
            task_set_intervals.emplace_back(t * sets_per_task, (t + 1) * sets_per_task - 1);
        }


        //create tasks to be handled by threads
        std::queue<DecompositionTask> decomp_tasks;
        for (const auto& t : task_set_intervals) {
            decomp_tasks.push((DecompositionTask){t, options.batch_size,data_name});
        }

        cout << "Starting Decomposition into connected sets using txt.\n\n";

        //start status thread
        std::atomic<int> processed_counter(0);
        std::thread status_thread(Threading::statusTHREAD, std::ref(processed_counter),num_components);

        //init mutexes
        std::mutex queue_mutex;
        std::condition_variable cv;
        std::atomic<bool> done_flag(false);

        //start worker threads
        std::vector<std::thread> decomp_threads;
        for(int t=0; t < num_threads; t++) {
            decomp_threads.emplace_back(Threading::decompositionFromTXTWorkerTHREAD,
                                        std::ref(decomp_tasks),std::ref(queue_mutex),std::ref(cv), std::ref(done_flag),
                                         std::ref(polys1), std::ref(polys2),std::ref(processed_counter));
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
        processed_counter = num_components;
        status_thread.join();

        cout << "Completed inital decomposition into connected subsets." << endl;
        std::chrono::steady_clock::time_point decomp_end = std::chrono::steady_clock::now();
        timing_decomp = std::chrono::duration_cast<std::chrono::nanoseconds>(decomp_end - decomp_begin).count() / 1e9;


    }
    else {
        throw std::runtime_error("No decomposition file found. Please provide a file with the decomposition of the input data as txt using 'matching_predecomposer.py'.");
    }

    //initialize one global Gurobi Environment to avoid collisions on license check in multiple threads
    GRBEnv env = GRBEnv(true);
    //deactivate console logging
    env.set("LogToConsole", "0");
    //start environment
    env.start();


    //execute for all lambdas
    for (const auto& lambda : options.lambdas) {
        //initialize logger
        Logger logger("../export/" + data_name + "/" + options.log_name + ".csv");
        logger.setInstanceName(data_name);
        logger.setLambda(lambda);
        logger.setExploitMode(options.exploit_opt_props);
        if (options.tree_mode == TREE_BUILD_MODE::INFORMED) logger.setTreeMode("informed");
        else logger.setTreeMode("kruskal");
        if (options.mode == SOLUTION_MODE::OPT) logger.setSolMode("opt");
        else if (options.mode == SOLUTION_MODE::CANZAR3APPROX) logger.setSolMode("canzar3approx");
        if (options.objective == OBJECTIVE::JACCARD) logger.setObjectiveMode("jaccard");
        else if (options.objective == OBJECTIVE::JACCARD_HAUSDORFF)
            logger.setObjectiveMode("jaccard_hausdorff_"+to_string(options.objective_weights.first) + "_" + to_string(options.objective_weights.second));
        logger.setInputPolygons(polys1.size(),polys2.size());
        logger.setConnectedComponents(merged_polys.size());



        //make sure the target folder exists
        //define subfolder export name depending on objective
        std::string target_subfolder_obj = "lambda" + to_string((int)(lambda * 100));
        if (options.objective == OBJECTIVE::JACCARD_HAUSDORFF) {
            target_subfolder_obj = "hd-jac-" + to_string((int)(options.objective_weights.first * 100))
                                    + "-" + to_string((int)(options.objective_weights.second * 100));
        }


        std::string target_folder = "../export/" + data_name + "/" + target_subfolder_obj;
        //if the jaccard hausdorff objective is chosen, change output folder

        if(!fs::exists(target_folder)) fs::create_directories(target_folder);

        //measure time of global execution
        std::chrono::steady_clock::time_point exec_begin = std::chrono::steady_clock::now();

        //prepare thread task pool
        std::queue<ConnectedSetTask> solving_tasks;
        for (int i=0; i< num_batches; i++) {
            solving_tasks.push((ConnectedSetTask){data_name, lambda, options,i});
        }


        cout << "Starting to solve for lambda: " << lambda << ".\n\n";
        //measure time
        std::chrono::steady_clock::time_point solving_begin = std::chrono::steady_clock::now();

        //start status thread
        std::atomic<int> processed_counter(0);
        auto status_thread = std::thread(Threading::statusTHREAD, std::ref(processed_counter), num_components);

        //init mutexes
        std::mutex queue_mutex;
        std::condition_variable cv;

        //start worker threads
        std::atomic<bool> done_flag(false);
        std::vector<std::thread> set_threads(num_threads);
        for(int t=0; t < num_threads; t++) {
            set_threads[t] = std::thread(Threading::solvingWorkerTHREAD,
                                         std::ref(solving_tasks),std::ref(queue_mutex),std::ref(cv),std::ref(done_flag),
                                         std::ref(env), std::ref(processed_counter));
        }

        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            done_flag = true;
        }
        cv.notify_all();

        //join threads
        for(int t=0; t < num_threads; t++) {
            set_threads[t].join();
        }

        // Ensure the status thread finishes
        processed_counter = num_components;
        status_thread.join();

        cout << "done, now combining solutions" << endl;

        //combine all part-solutions produced by the threads
        std::vector<int> set_ids;
        std::vector<std::pair<int,int>> set_sizes;
        std::vector<std::vector<double>> execution_times;
        Solution sol = Threading::combineThreadSolutions(data_name,num_batches,polys1.size(),polys2.size(),set_ids,
                                set_sizes,execution_times);

        cout << "completing matching..." << endl;

        //complete global matching (setting negative match IDs for unmatched polygons)
        sol.completeMatching();

        //this is where the pure computation ends, take time and write to logger
        std::chrono::steady_clock::time_point exec_end = std::chrono::steady_clock::now();
        cout << "Logging..." << endl;
        std::vector<double> timings; timings.resize(6,0.0);
        timings[0] = timing_decomp;
        timings[5] = std::chrono::duration_cast<std::chrono::nanoseconds>(exec_end - exec_begin).count() / 1e9;
        for (int i=1; i<5; i++) {
            for (int j=0; j< execution_times.size(); j++) {
                timings[i] += execution_times[j][i-1];
            }
        }
        logger.setTimings(timings);
        logger.setMemoryUsage(get_peak_memory_kb());
        logger.setObjective(sol.getTargetValue());
        logger.setTotalNumMatches(sol.getMatchCount());
        logger.setMatchDistribution(sol.getMatchDistribution());
        logger.setObjectiveDistribution(sol.getObjectiveDistribution());
        logger.log();
        std::cout << "\x1b[1A\x1b[2K";
        cout << "Logging completed." << endl;

        //export the global solution
        std::string output_path1 = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "_" + data_descr1 + "_matched";
        std::string output_path2 = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "_" + data_descr2 + "_matched";
        std::string csv_output_path = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "_data";
        std::string analysis_output_path = "../export/" + data_name + "/" + target_subfolder_obj + "/" + data_name + "_analysis";

        //if geopackage input, copy input files and extend attributes, else write shapefiles
        for (int i=0; i<2; i++) {
            std::string output_path = i==0? output_path1 + ".gpkg" : output_path2 + ".gpkg";
            fs::path gpkgPath = i==0? file1 + ".gpkg" : file2 + ".gpkg";
            const std::vector<Polygon_wh>& polys = i==0? polys1 : polys2;
            if (fs::exists(gpkgPath)) {
                std::filesystem::copy(gpkgPath, output_path, std::filesystem::copy_options::overwrite_existing);
                updateGeoPackageWithAttributes(output_path,polys,sol,i);
            } else {
                //write labeled polygons to shapefiles
                throw std::runtime_error("No geopackage file found, please provide geopackages as input.");
            }
        }


        //write analysis data to csvs
        writeToCSV(sol, csv_output_path);


        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "Computed n:m matching and exported results. (" << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [s])" << endl;
        cout << "Objective = " << sol.getTargetValue() << endl;
    }

    //cleanup tmp storage
    Threading::cleanup_temp_folder("../input/" + data_name + "/tmp/");

}

int main(int argc, char* argv[]) {
    try {
        CommandLineOptions options = parse_command_line(argc,argv);

        if (options.dataset_name.empty() ||options.lambdas.empty()) {
            std::cerr << "Error: Both -d and -l (or -lr) options are required." << endl;
            print_help();
            return 1;
        }

        run_NtoM_with_PreDecomposition(options);

        cout << "Completed matching " << options.dataset_name << " for lambda(s) ";
        for (const auto& l : options.lambdas) {cout << l << ", ";} cout << endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << endl;
        print_help();
        return 1;
    }


	return 0;
}
