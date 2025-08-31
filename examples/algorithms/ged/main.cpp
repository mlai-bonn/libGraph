// This file shows how to use the libGraph graph library


#include "SimplePatterns.h"
#include "Algorithms/GED/GEDFunctions.h"
#include "GraphDataStructures/GraphBase.h"
#include "src/env/ged_env.hpp"

// main function to run the examples

void add_graph_to_env(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, const GraphStruct& g) {
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

GEDResult result_from_env(ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID>& env, GraphData<GraphStruct>& graphs, const int source_graph_id = 0, const int target_graph_id = 1){
    const ged::NodeMap node_map = env.get_node_map(source_graph_id, target_graph_id);
    std::pair<Nodes, Nodes> mapping;
    for (const auto x : node_map.get_forward_map()) {
        mapping.first.push_back(x);
    }
    for (const auto x : node_map.get_backward_map()) {
        mapping.second.push_back(x);
    }

    GEDResult result = {
        node_map.induced_cost(),
        {graphs[source_graph_id], graphs[target_graph_id]},
        mapping,
        env.get_runtime(source_graph_id,target_graph_id),
    };
    return result;
}

int main() {

    const GraphStruct triangle = SimplePatterns::Circle(3);
    const GraphStruct square = SimplePatterns::Circle(5);

    GraphData<GraphStruct> graph_data;
    graph_data.add(triangle);
    graph_data.add(square);

    ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID> env;
    env.set_edit_costs(ged::Options::EditCosts::CONSTANT);
    add_graph_to_env(env, triangle);
    add_graph_to_env(env, square);
    env.init();


    env.set_method(ged::Options::GEDMethod::STAR);
    env.init_method();
    constexpr int source_id = 1;
    constexpr int target_id = 0;
    env.run_method(source_id,target_id);


    GEDResult result = result_from_env(env, graph_data, source_id, target_id);
    std::cout << "Approximated Distance: " << result.distance << std::endl;
        std::cout << "Time: " << result.time << " seconds" << std::endl;
    std::cout << "Quasimetric Cost: " << env.quasimetric_costs() << std::endl;
    EditPath edit_path;
    result.get_edit_path(edit_path, 0);
    GraphData<GraphStruct> edit_path_graphs;
    for (const auto& g : edit_path.edit_path_graphs) {
        edit_path_graphs.add(g);
    }
    edit_path_graphs.graphData.back().SetName(edit_path.target_graph.GetName());
    // print all the graphs
    for (const auto& g : edit_path_graphs.graphData) {
        std::cout << g << std::endl;
    }
    return 0;


//    algorithms::GEDDFS solver(source, target, nullptr, nullptr, 0, 0, false, false);

//    test();
//    GEDApproximationParameters parameters = graph_edit_distance_approximation();
//    return 0;
}

