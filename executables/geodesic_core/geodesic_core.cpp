


#include <string>
#include <filesystem>
#include <vector>
#include "../../include/Algorithms/GeodesicCore/CoreGrowAlgorithm.h"
#include "../../include/DataClasses.h"
#include "../../include/Algorithms/GeodesicCore/CoreAlgorithm.h"
#include "../../include/io/FileEvaluation.h"

void grow_vs_intersection_core()
{
    // num runs
    int num_runs = 100;
    int num_steps = 15;
    double threshold = 0.9;
    // iterate over several graphs
    std::string directory = "Collaboration";
    // paths
    std::string directory_path = "../../../../../GraphData/RealWorld/" + directory + "/";
    // get all files with extension .edges in the directory using std::filesystem
    std::vector<std::string> graph_paths;
    for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
        if (entry.path().extension() == ".edges") {
            // check if entry path contains "CA-GrQc"
            if (entry.path().string().find("CA-GrQc") == std::string::npos) {
                continue;
            }
            graph_paths.emplace_back(entry.path().string());
        }
    }

    for (const auto &path: graph_paths) {

        // load graph
        GraphStruct graph = GraphStruct(path);
        // print the graph information (nodes, edges, density)
        std::cout << "Nodes: " << graph.nodes() << std::endl;
        std::cout << "Edges: " << graph.edges() << std::endl;
        std::cout << "Density: " << (double) graph.edges()/(graph.nodes()*graph.nodes())  << std::endl;


        // run the grow core algorithm
        auto core_grow_input_paramaters = CoreGrowAlgorithmInputParameters{num_runs, num_steps, threshold, 0,
                                                                           true, true};
        CoreGrowAlgorithmOutputParameters core_grow_output_parameters;
        // run the core grow algorithm
        // measure the runtime
        auto start = std::chrono::high_resolution_clock::now();
        CoreGrowAlgorithm CGA = CoreGrowAlgorithm(graph);
        CGA.Run(core_grow_output_parameters, core_grow_input_paramaters);
        // print the core size
        std::cout << "Grow Core size: " << core_grow_output_parameters.core_nodes.size() << std::endl;

        // run our core algorithm
        auto core_input_paramaters = CoreAlgorithmInputParameters{.generator_size = num_steps, .print = true, .save = true};
        CoreAlgorithmOutputParameters core_output_parameters;
        CoreAlgorithm coreAlgorithm = CoreAlgorithm(graph);
        coreAlgorithm.Run(core_output_parameters, core_input_paramaters);
        // print the core size
        std::cout << "Core size: " << core_output_parameters.core_nodes.size() << std::endl;
        
        // compute Jaccard similarity between core_grow_output_parameters.core_nodes and coreNodes i.e. the size of the intersection divided by the size of the union
        // intersection of core_grow_output_parameters.core_nodes and coreNodes
        std::vector<NodeId> _intersection;
        // union of core_grow_output_parameters.core_nodes and coreNodes
        std::vector<NodeId> _union;
        // compute the intersection and the union
        std::set_intersection(core_grow_output_parameters.core_nodes.begin(), core_grow_output_parameters.core_nodes.end(),
                              core_output_parameters.core_nodes.begin(), core_output_parameters.core_nodes.end(), std::back_inserter(_intersection));
        std::set_union(core_grow_output_parameters.core_nodes.begin(), core_grow_output_parameters.core_nodes.end(), core_output_parameters.core_nodes.begin(),
                       core_output_parameters.core_nodes.end(), std::back_inserter(_union));
        // compute the Jaccard similarity
        double jaccard_similarity = (double) _intersection.size() / (double) _union.size();
        // print the Jaccard similarity
        std::cout << "Jaccard similarity: " << jaccard_similarity << std::endl;

        // save results to file using FileEvaluation
        FileEvaluation fileEvaluation = FileEvaluation("../results/", directory, ".csv");
        fileEvaluation.headerValueInsert(
                {"Graph", "Nodes", "Edges", "Density", "Grow Core Size", "Grow Runtime", "Exact Core Size", "Exact Runtime", "Jaccard similarity"},
                {graph.GetName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                 std::to_string((double) graph.edges()/(graph.nodes()*graph.nodes())), std::to_string(core_grow_output_parameters.core_nodes.size()), std::to_string(core_grow_output_parameters.runtime),
                 std::to_string(core_output_parameters.core_nodes.size()), std::to_string(core_output_parameters.runtime), std::to_string(jaccard_similarity)});
        fileEvaluation.save();
    }

}

void stability_grow_core() {
    // num runs
    auto num_runs_list = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    auto num_steps_list = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};

    for (auto num_runs: num_runs_list) {
        for (auto num_steps: num_steps_list) {

            // iterate over several graphs
            std::string directory = "Collaboration";
            // paths
            std::string directory_path = "../../../../../GraphData/RealWorld/" + directory + "/";
            // get all files with extension .edges in the directory using std::filesystem
            std::vector<std::string> graph_paths;
            for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
                if (entry.path().extension() == ".edges") {
                    // check if entry path contains "CA-GrQc"
                    if (entry.path().string().find("CA-GrQc") == std::string::npos) {
                        continue;
                    }
                    graph_paths.emplace_back(entry.path().string());
                }
            }

            for (const auto &path: graph_paths) {
                auto inputParameters = CoreGrowAlgorithmInputParameters{num_runs, num_steps, 0.9, 0,
                                                                        true, true};
                CoreGrowAlgorithmOutputParameters outputParameters;
                GraphStruct graph = GraphStruct(path);
                // print the graph information (nodes, edges, density)
                std::cout << "Nodes: " << graph.nodes() << std::endl;
                std::cout << "Edges: " << graph.edges() << std::endl;
                std::cout << "Density: " << (double) graph.edges()/(graph.nodes()*graph.nodes()) << std::endl;
                // run the core grow algorithm
                // measure the runtime
                auto start = std::chrono::high_resolution_clock::now();
                auto CGA = CoreGrowAlgorithm(graph);
                //CGA.Run(outputParameters, inputParameters);
                double grow_runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - start).count() /
                                       1000000.0);
                // print the core size
                std::cout << "Core size: " << outputParameters.core_nodes.size() << std::endl;

                // save results to file using FileEvaluation
                FileEvaluation fileEvaluation = FileEvaluation("../results/", directory + "_stability_grow", ".csv");
                fileEvaluation.headerValueInsert(
                        {"Graph", "Nodes", "Edges", "Density", "Grow Core Size", "Grow Runtime", "Num Runs",
                         "Num Steps"},
                        {graph.GetName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                         std::to_string((double) graph.edges()/(graph.nodes()*graph.nodes())), std::to_string(outputParameters.core_nodes.size()),
                         std::to_string(grow_runtime), std::to_string(num_runs), std::to_string(num_steps)});
                fileEvaluation.save();
            }
        }
    }
}

void stabilility_intersection_core() {
    // num runs
    auto num_steps_list = {5};//, 10, 15, 20, 25, 30, 35, 40, 45, 50};


    // iterate over several graphs
    std::string directory = "Collaboration";
    // paths
    std::string directory_path = "../../../../../GraphData/RealWorld/" + directory + "/";
    // get all files with extension .edges in the directory using std::filesystem
    std::vector<std::string> graph_paths;
    for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
        if (entry.path().extension() == ".edges") {
            // check if entry path contains "CA-GrQc"
            if (entry.path().string().find("CA-GrQc") == std::string::npos) {
                continue;
            }
            graph_paths.emplace_back(entry.path().string());
        }
    }

    for (const auto &path: graph_paths) {
        GraphStruct graph = GraphStruct(path);
        // print the graph information (nodes, edges, density)
        std::cout << "Nodes: " << graph.nodes() << std::endl;
        std::cout << "Edges: " << graph.edges() << std::endl;
        std::cout << "Density: " << (double) graph.edges() / (graph.nodes() * graph.nodes()) << std::endl;

        for (auto num_steps: num_steps_list) {
            // run our core algorithm
            auto core_input_paramaters = CoreAlgorithmInputParameters{.generator_size = num_steps, .print = true, .save = true};
            CoreAlgorithmOutputParameters core_output_parameters;
            CoreAlgorithm coreAlgorithm = CoreAlgorithm(graph);
            coreAlgorithm.Run(core_output_parameters, core_input_paramaters);
            // print the core size
            std::cout << "Core size: " << core_output_parameters.core_nodes.size() << std::endl;


            // save results to file using FileEvaluation
            FileEvaluation fileEvaluation = FileEvaluation("../results/", directory + "_stability_intersection",
                                                           ".csv");
            fileEvaluation.headerValueInsert(
                    {"Graph", "Nodes", "Edges", "Density", "Intersection Core Size", "Exact Runtime", "Num Steps"},
                    {graph.GetName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                     std::to_string((double) graph.edges() / (graph.nodes() * graph.nodes())),
                     std::to_string(core_output_parameters.core_nodes.size()),
                     std::to_string(core_output_parameters.runtime), std::to_string(num_steps)});
            fileEvaluation.save();
        }
    }
}


int main(int argc, char *argv[]) {
    grow_vs_intersection_core();
    //stability_grow_core();
    //stabilility_intersection_core();
}
