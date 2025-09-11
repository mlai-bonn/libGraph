//
// Created by florian on 09.09.25.
//

#ifndef GEDEXAMPLE_GEDLIBWRAPPER_H
#define GEDEXAMPLE_GEDLIBWRAPPER_H

// For the following code GEDLIB has to be included and #define GEDLIB must be set in the main file
#ifdef GEDLIB
#include "GEDFunctions.h"
#include "src/env/ged_env.hpp"

// description of the wrapper functions

/* Adds a graph to the given GEDLIB environment
 *
 * @param env The GEDLIB environment
 * @param g The graph to be added
 */
template<typename T>
void AddGraphToGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const T& g);
/* Adds multiple graphs to the given GEDLIB environment
 * @param env The GEDLIB environment
 * @param graph_data The graphs to be added
 */
template<typename T>
void AddGraphsToGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graph_data);
/* Initializes the given GEDLIB environment
 * @param env The GEDLIB environment
 * @param graph_data The graphs to be added
 * @param edit_costs The edit costs to be used
 * @param method The method to be used
 */
template<typename T>
void InitializeGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graph_data, ged::Options::EditCosts edit_costs = ged::Options::EditCosts::CONSTANT, ged::Options::GEDMethod method = ged::Options::GEDMethod::BP_BEAM);
/* Evaluates the result of a GEDLIB computation
 * @param env The GEDLIB environment
 * @param graphs The graphs used in the computation
 * @param source_graph_id The ID of the source graph
 * @param target_graph_id The ID of the target graph
 * @return The evaluation of the result
 */
template<typename T>
GEDEvaluation<T> ComputeGEDResult(const ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, GraphData<T>& graphs, const int source_graph_id = 0, const int target_graph_id = 1);
/* Evaluates all results of a GEDLIB computation and saves the mappings to a file
 * @param env The GEDLIB environment
 * @param graph_data The graphs used in the computation
 * @param mapping_path_output The path to the file where the mappings should be saved
 */
template<typename T>
void ComputeGEDResults(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID> &env,
                                   const GraphData<T> &graph_data,
                                   const std::vector<std::pair<INDEX,INDEX>> &graph_ids,
                                   const std::string &results_path);



template<typename T>
void AddGraphToGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const T& g) {
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
        for (const auto j : g.get_neighbors(i)) {
            if (i < j) {
                env.add_edge(graph_ids.second - 1, i, j,0);
            }
        }
    }
}

template<typename T>
void AddGraphsToGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graph_data) {
    for (const auto& g : graph_data.graphData) {
        AddGraphToGEDEnvironment(env, g);
    }
}

template<typename T>
void InitializeGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graph_data, ged::Options::EditCosts edit_costs, ged::Options::GEDMethod method) {
    env.set_edit_costs(edit_costs);
    AddGraphsToGEDEnvironment(env, graph_data);
    env.init();
    env.set_method(method);
    env.init_method();
}

template<typename T>
GEDEvaluation<T> ComputeGEDResult(const ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graphs, const int source_graph_id, const int target_graph_id){
    ged::NodeMap node_map = env.get_node_map(source_graph_id, target_graph_id);
    std::pair<Nodes, Nodes> mapping;
    for (const auto x : node_map.get_forward_map()) {
        mapping.first.push_back(x);
    }
    for (const auto x : node_map.get_backward_map()) {
        mapping.second.push_back(x);
    }
    ged::GEDGraph::GraphID i = source_graph_id;
    ged::GEDGraph::GraphID j = target_graph_id;
    env.compute_induced_cost(i, j, node_map);
    GEDEvaluation<T> result = {
        env.get_node_map(i,j).induced_cost(),
        env.get_lower_bound(i,j),
        env.get_upper_bound(i,j),
        {graphs.graphData[source_graph_id], graphs.graphData[target_graph_id]},
        {source_graph_id, target_graph_id},
        mapping,
        graphs.GetName(),
        env.get_runtime(source_graph_id,target_graph_id),
    };
    return result;
}

template<typename T>
void ComputeGEDResults(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID> &env,
                                   const GraphData<T> &graph_data,
                                   const std::vector<std::pair<INDEX,INDEX>> &graph_ids,
                                   const std::string &results_path) {
    // check whether the output path exists
    if (!std::filesystem::exists(results_path)) {
        std::cout << "The output path " << results_path << " does not exist." << std::endl;
        return;
    }
    // check whether the file already exists in the output_path
    if (std::filesystem::exists(results_path + graph_data.GetName() + "_ged_mapping.bin")) {
        std::cout << "The mapping file for " << graph_data.GetName() << " already exist." << std::endl;
        std::cout << "Skipping computation." << std::endl;
        return;
    }
    // create tmp directory
    std::filesystem::create_directory(results_path + "tmp/");
    // counter for number of computed results
    int counter = 0;
    // time variable
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (const auto& [i,j] : graph_ids) {
        if (j > i) {
            // Check if mapping already exists in the tmp folder
            if (std::filesystem::exists(results_path + "tmp/" + graph_data.GetName() + "_" + std::to_string(i) + "_" + std::to_string(j) + "_ged_mapping.bin")) {
                //std::cout << "Mapping between graph " << i << " and graph " << j << " already exists. Skipping." << std::endl;
                ++counter;
                continue;
            }
            //std::cout << "Computing mapping between graph " << i << " and graph " << j << std::endl;
            // print percentage
            //std::cout << "Progress: " << (counter * 100) / (graph_data.size() * (graph_data.size() - 1) / 2) << "%" << std::endl;
            // estimated time in minutes
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            const double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            const double estimated_total_time = (elapsed_seconds / (counter + 1)) * (graph_data.size() * (graph_data.size() - 1) / 2);
            const double estimated_time_left = estimated_total_time - elapsed_seconds;
            //std::cout << "Estimated time left: " << estimated_time_left / 60 << " minutes" << std::endl;
            env.run_method(i, j);
            GEDEvaluation<T> result = ComputeGEDResult(env, graph_data, i, j);
            // print calculated (approximated Distance)
            //std::cout << "Approximated Distance " << i << " to " << j << " : " << result.distance << std::endl;
            // save result to binary
            GEDResultToBinary(results_path + "tmp/", result);
            //std::cout << "Saved intermediate result for graphs " << i << " and " << j << std::endl;
            ++counter;
        }
    }
}





#endif

#endif //GEDEXAMPLE_GEDLIBWRAPPER_H