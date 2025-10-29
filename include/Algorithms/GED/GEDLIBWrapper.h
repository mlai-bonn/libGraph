//
// Created by florian on 09.09.25.
//

#ifndef GED_LIB_WRAPPER_H
#define GED_LIB_WRAPPER_H

// For the following code GEDLIB has to be included and #define GEDLIB must be set in the main file
#ifdef GEDLIB
#include "GEDFunctions.h"
#include "GEDEvaluation.h"
#include "src/env/ged_env.hpp"
#include <numeric>

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
void InitializeGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graph_data, ged::Options::EditCosts edit_costs = ged::Options::EditCosts::CONSTANT, ged::Options::GEDMethod method = ged::Options::GEDMethod::BP_BEAM, const std::string& method_options = "");
/* Evaluates the result of a GEDLIB computation
 * @param env The GEDLIB environment
 * @param graphs The graphs used in the computation
 * @param source_graph_id The ID of the source graph
 * @param target_graph_id The ID of the target graph
 * @return The evaluation of the result
 */
template<typename T>
GEDEvaluation<T> ComputeGEDResult(const ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graphs, const int source_graph_id = 0, const int target_graph_id = 1);
/* Evaluates all results of a GEDLIB computation and saves the mappings to a file
 * @param env The GEDLIB environment
 * @param graph_data The graphs used in the computation
 * @param mapping_path_output The path to the file where the mappings should be saved
 */
template<typename T>
void ComputeGEDResults(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID> &env,
                                   const GraphData<T> &graph_data,
                                   const std::vector<std::pair<INDEX,INDEX>> &graph_ids,
                                   const std::string &results_path, ged::Options::GEDMethod method, const std::string& method_options);





template<typename T>
void AddGraphToGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const T& g) {
    // T must be from base class GraphStruct
    static_assert(std::is_base_of<GraphStruct, T>::value, "T must be from base class GraphStruct");
    // Add graph
    env.add_graph(g.GetName());
    // Add nodes
    // get last env graph id
    ged::GEDGraph::GraphID last_graph_id = env.graph_ids().second - 1;
    for (int i = 0; i < g.nodes(); ++i) {
        if (g.labelType == LABEL_TYPE::UNLABELED) {
            env.add_node(last_graph_id, i, 0);
        }
        else {
            env.add_node(last_graph_id, i, g.label(i));
        }
    }
    // Add edges
    for (int i = 0; i < g.nodes(); ++i) {
        for (const auto j : g.get_neighbors(i)) {
            if (i < j) {
                INDEX edge_label = (INDEX) g.GetEdgeData({i,j}, "label");
                env.add_edge(last_graph_id, i, j, edge_label);
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
void InitializeGEDEnvironment(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphData<T>& graph_data, ged::Options::EditCosts edit_costs, ged::Options::GEDMethod method, const std::string& method_options) {
    env.set_edit_costs(edit_costs);
    AddGraphsToGEDEnvironment(env, graph_data);
    env.init();
    env.set_method(method, method_options);
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
                                   const std::string &results_path, ged::Options::GEDMethod method, const std::string& method_options) {
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
    std::filesystem::create_directory(results_path);
    // counter for number of computed results
    int counter = 0;
    // time variable
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (const auto& [i,j] : graph_ids) {
        if (j > i) {
            // Check if mapping already exists in the tmp folder
            if (std::filesystem::exists(results_path + graph_data.GetName() + "_" + std::to_string(i) + "_" + std::to_string(j) + "_ged_mapping.bin")) {
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
            // check if the mapping is valid
            bool valid = CheckResultsValidity(std::vector<GEDEvaluation<T>>{result}).empty();
            // Try to fix using only one thread
            if (!valid) {
                std::cout << "The computation of the mapping between graph " << i << " and graph " << j << " resulted in an invalid mapping. Retrying with single thread." << std::endl;
                // set --threads 1 in method
                ged::Options::GEDMethod old_method = method;
                // replace --threads x with --threads 1
                std::string new_method_options = method_options;
                const size_t pos = new_method_options.find("--threads");;
                if (pos != std::string::npos) {
                    size_t end_pos = new_method_options.find(" ", pos + 9);
                    if (end_pos == std::string::npos) {
                        end_pos = new_method_options.length();
                    }
                    new_method_options.replace(pos, end_pos - pos, "--threads 1");
                }
                else {
                    if (!new_method_options.empty() && new_method_options.back() != ' ') {
                        new_method_options += " ";
                    }
                    new_method_options += "--threads 1";
                }
                env.set_method(old_method, new_method_options);
                env.run_method(i, j);
                result = ComputeGEDResult(env, graph_data, i, j);
                valid = CheckResultsValidity(std::vector<GEDEvaluation<T>>{result}).empty();
                if (valid) {
                    std::cout << "Successfully fixed the mapping between graph " << i << " and graph " << j << " by using a single thread." << std::endl;
                }
                else {
                    std::cout << "Failed to fix the mapping between graph " << i << " and graph " << j << "." << std::endl;
                }
                // reset method
                env.set_method(method, method_options);
            }
            // print calculated (approximated Distance)
            //std::cout << "Approximated Distance " << i << " to " << j << " : " << result.distance << std::endl;
            // save result to binary
            GEDResultToBinary(results_path, result);
            //std::cout << "Saved intermediate result for graphs " << i << " and " << j << std::endl;
            ++counter;
        }
    }
}

std::string GEDMethodToString(const ged::Options::GEDMethod method) {
    switch (method) {
#ifdef GUROBI
        case ged::Options::GEDMethod::F1:
            return "F1";
        case ged::Options::GEDMethod::F2:
            return "F2";
        case ged::Options::GEDMethod::COMPACT_MIP:
            return "COMPACT_MIP";
        case ged::Options::GEDMethod::BLP_NO_EDGE_LABELS:
            return "BLP_NO_EDGE_LABELS";
#endif
        case ged::Options::GEDMethod::BRANCH:
            return "BRANCH";
        case ged::Options::GEDMethod::BRANCH_FAST:
            return "BRANCH_FAST";
        case ged::Options::GEDMethod::BRANCH_TIGHT:
            return "BRANCH_TIGHT";
        case ged::Options::GEDMethod::BRANCH_UNIFORM:
            return "BRANCH_UNIFORM";
        case ged::Options::GEDMethod::BRANCH_COMPACT:
            return "BRANCH_COMPACT";
        case ged::Options::GEDMethod::PARTITION:
            return "PARTITION";
        case ged::Options::GEDMethod::HYBRID:
            return "HYBRID";
        case ged::Options::GEDMethod::RING:
            return "RING";
        case ged::Options::GEDMethod::ANCHOR_AWARE_GED:
            return "ANCHOR_AWARE_GED";
        case ged::Options::GEDMethod::WALKS:
            return "WALKS";
        case ged::Options::GEDMethod::IPFP:
            return "IPFP";
        case ged::Options::GEDMethod::BIPARTITE:
            return "BIPARTITE";
        case ged::Options::GEDMethod::SUBGRAPH:
            return "SUBGRAPH";
        case ged::Options::GEDMethod::NODE:
            return "NODE";
        case ged::Options::GEDMethod::RING_ML:
            return "RING_ML";
        case ged::Options::GEDMethod::BIPARTITE_ML:
            return "BIPARTITE_ML";
        case ged::Options::GEDMethod::REFINE:
            return "REFINE";
        case ged::Options::GEDMethod::BP_BEAM:
            return "BP_BEAM";
        case ged::Options::GEDMethod::SIMULATED_ANNEALING:
            return "SIMULATED_ANNEALING";
        case ged::Options::GEDMethod::HED:
            return "HED";
        case ged::Options::GEDMethod::STAR:
            return "STAR";
        default:
            return "Unknown";
    }
}

ged::Options::GEDMethod GEDMethodFromString(const std::string& methodStr) {
#ifdef GUROBI
    if (methodStr == "F1") return ged::Options::GEDMethod::F1;
    if (methodStr == "F2") return ged::Options::GEDMethod::F2;
    if (methodStr == "COMPACT_MIP") return ged::Options::GEDMethod::COMPACT_MIP;
    if (methodStr == "BLP_NO_EDGE_LABELS") return ged::Options::GEDMethod::BLP_NO_EDGE_LABELS;
#endif
    if (methodStr == "BRANCH") return ged::Options::GEDMethod::BRANCH;
    if (methodStr == "BRANCH_FAST") return ged::Options::GEDMethod::BRANCH_FAST;
    if (methodStr == "BRANCH_TIGHT") return ged::Options::GEDMethod::BRANCH_TIGHT;
    if (methodStr == "BRANCH_UNIFORM") return ged::Options::GEDMethod::BRANCH_UNIFORM;
    if (methodStr == "BRANCH_COMPACT") return ged::Options::GEDMethod::BRANCH_COMPACT;
    if (methodStr == "PARTITION") return ged::Options::GEDMethod::PARTITION;
    if (methodStr == "HYBRID") return ged::Options::GEDMethod::HYBRID;
    if (methodStr == "RING") return ged::Options::GEDMethod::RING;
    if (methodStr == "ANCHOR_AWARE_GED") return ged::Options::GEDMethod::ANCHOR_AWARE_GED;
    if (methodStr == "WALKS") return ged::Options::GEDMethod::WALKS;
    if (methodStr == "IPFP") return ged::Options::GEDMethod::IPFP;
    if (methodStr == "BIPARTITE") return ged::Options::GEDMethod::BIPARTITE;
    if (methodStr == "SUBGRAPH") return ged::Options::GEDMethod::SUBGRAPH;
    if (methodStr == "NODE") return ged::Options::GEDMethod::NODE;
    if (methodStr == "RING_ML") return ged::Options::GEDMethod::RING_ML;
    if (methodStr == "BIPARTITE_ML") return ged::Options::GEDMethod::BIPARTITE_ML;
    if (methodStr == "REFINE") return ged::Options::GEDMethod::REFINE;
    if (methodStr == "BP_BEAM") return ged::Options::GEDMethod::BP_BEAM;
    if (methodStr == "SIMULATED_ANNEALING") return ged::Options::GEDMethod::SIMULATED_ANNEALING;
    if (methodStr == "HED") return ged::Options::GEDMethod::HED;
    if (methodStr == "STAR") return ged::Options::GEDMethod::STAR;
    // Fallback: BP_BEAM
    return ged::Options::GEDMethod::BP_BEAM;
}

std::string EditCostsToString(const ged::Options::EditCosts costs) {
    switch (costs) {
        case ged::Options::EditCosts::CHEM_1:
            return "CHEM_1";
        case ged::Options::EditCosts::CHEM_2:
            return "CHEM_2";
        case ged::Options::EditCosts::CMU:
            return "CMU";
        case ged::Options::EditCosts::GREC_1:
            return "GREC_1";
        case ged::Options::EditCosts::GREC_2:
            return "GREC_2";
        case ged::Options::EditCosts::PROTEIN:
            return "PROTEIN";
        case ged::Options::EditCosts::FINGERPRINT:
            return "FINGERPRINT";
        case ged::Options::EditCosts::LETTER:
            return "LETTER";
        case ged::Options::EditCosts::CONSTANT:
            return "CONSTANT";
        default:
            return "Unknown";
    }
}

ged::Options::EditCosts EditCostsFromString(const std::string& costsStr) {
    if (costsStr == "CHEM_1") return ged::Options::EditCosts::CHEM_1;
    if (costsStr == "CHEM_2") return ged::Options::EditCosts::CHEM_2;
    if (costsStr == "CMU") return ged::Options::EditCosts::CMU;
    if (costsStr == "GREC_1") return ged::Options::EditCosts::GREC_1;
    if (costsStr == "GREC_2") return ged::Options::EditCosts::GREC_2;
    if (costsStr == "PROTEIN") return ged::Options::EditCosts::PROTEIN;
    if (costsStr == "FINGERPRINT") return ged::Options::EditCosts::FINGERPRINT;
    if (costsStr == "LETTER") return ged::Options::EditCosts::LETTER;
    if (costsStr == "CONSTANT") return ged::Options::EditCosts::CONSTANT;
    // Fallback: CONSTANT
    return ged::Options::EditCosts::CONSTANT;
}





#endif

#endif //GED_LIB_WRAPPER_H
