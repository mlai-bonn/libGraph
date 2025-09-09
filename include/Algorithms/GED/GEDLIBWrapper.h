//
// Created by florian on 09.09.25.
//

#ifndef GEDEXAMPLE_GEDLIBWRAPPER_H
#define GEDEXAMPLE_GEDLIBWRAPPER_H

// For the following code GEDLIB has to be included and #define GEDLIB must be set in the main file
#ifdef GEDLIB
#include "GEDFunctions.h"
#include "src/env/ged_env.hpp"

inline void AddGraphToGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphStruct& g) {
    // Add graph
    env.add_graph(g.GetName());
    // Add nodes
    // get last env graph id
    std::pair<ged::GEDGraph::GraphID, ged::GEDGraph::GraphID> graph_ids = env.graph_ids();
    for (int i = 0; i < g.nodes(); ++i) {
        if (g.labelType == LABEL_TYPE::UNLABELED) {
            env.add_node(graph_ids.second - 1, i, 0);
        }
        else {
            env.add_node(graph_ids.second - 1, i, g.label(i));
        }
    }
    // Add edges
    for (int i = 0; i < g.nodes(); ++i) {
        for (int j = i + 1; j < g.nodes(); ++j) {
            env.add_edge(graph_ids.second - 1, i, j,0);
        }
    }
}

inline void AddGraphsToGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<GraphStruct>& graph_data) {
    for (const auto& g : graph_data.graphData) {
        AddGraphToGEDEnvironment(env, g);
    }
}

inline void InitializeGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<GraphStruct>& graph_data) {
   AddGraphsToGEDEnvironment(env, graph_data);
    env.init();
}

inline GEDEvaluation EvaluateGEDResult(const ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, GraphData<GraphStruct>& graphs, const int source_graph_id = 0, const int target_graph_id = 1){
    const ged::NodeMap& node_map = env.get_node_map(source_graph_id, target_graph_id);
    std::pair<Nodes, Nodes> mapping;
    for (const auto x : node_map.get_forward_map()) {
        mapping.first.push_back(x);
    }
    for (const auto x : node_map.get_backward_map()) {
        mapping.second.push_back(x);
    }

    GEDEvaluation result = {
        node_map.induced_cost(),
        {graphs[source_graph_id], graphs[target_graph_id]},
        {source_graph_id, target_graph_id},
        mapping,
        graphs.GetName(),
        env.get_runtime(source_graph_id,target_graph_id),
    };
    return result;
}

inline void EvaluateGEDResults(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID> &env,
                                   GraphData<GraphStruct> &graph_data,
                                   const std::string &mapping_path_output = "../Data/Mappings/") {
    // check whether file already exists
    if (std::filesystem::exists(mapping_path_output + graph_data.GetName() + "_ged_mapping.bin")) {
        std::cout << "The mapping file for " << graph_data.GetName() << " already exist." << std::endl;
        return;
    }
    // create tmp directory
    std::filesystem::create_directory(mapping_path_output + "tmp/");
    // counter for number of computed paths
    int counter = 0;
    // time variable
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (int i = 0; i < graph_data.size(); ++i) {
        for (int j = i + 1; j < graph_data.size(); ++j) {
            // Check if mapping already exists in the tmp folder
            if (std::filesystem::exists(mapping_path_output + "tmp/" + graph_data.GetName() + "_" + std::to_string(i) + "_" + std::to_string(j) + "_ged_mapping.bin")) {
                std::cout << "Mapping between graph " << i << " and graph " << j << " already exists. Skipping." << std::endl;
                ++counter;
                continue;
            }
            std::cout << "Computing Path between graph " << i << " and graph " << j << std::endl;
            // print percentage
            std::cout << "Progress: " << (counter * 100) / (graph_data.size() * (graph_data.size() - 1) / 2) << "%" << std::endl;
            // estimated time in minutes
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            const double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            const double estimated_total_time = (elapsed_seconds / (counter + 1)) * (graph_data.size() * (graph_data.size() - 1) / 2);
            const double estimated_time_left = estimated_total_time - elapsed_seconds;
            std::cout << "Estimated time left: " << estimated_time_left / 60 << " minutes" << std::endl;
            std::cout << "Computing Mapping between graph " << i << " and graph " << j << std::endl;
            env.run_method(i, j);
            GEDEvaluation result = EvaluateGEDResult(env, graph_data, i, j);
            // save result to binary
            GEDResultToBinary(mapping_path_output + "tmp/", result);
            std::cout << "Saved intermediate result for graphs " << i << " and " << j << std::endl;
            ++counter;
        }
    }
    // Merge all mappings in tmp folder
    std::vector<GEDEvaluation> merged_results;
    MergeBinaries(mapping_path_output + "tmp/", graph_data, merged_results);
    // Save final result
    GEDResultToBinary(mapping_path_output, merged_results);
    // remove all databasename related files in tmp folder
    for (const auto& entry : std::filesystem::directory_iterator(mapping_path_output + "tmp/")) {
        if (entry.path().string().find(graph_data.GetName() + "_") != std::string::npos) {
            std::filesystem::remove(entry.path());
        }
    }
}



#endif

#endif //GEDEXAMPLE_GEDLIBWRAPPER_H