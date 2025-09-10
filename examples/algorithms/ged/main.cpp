// This file shows how to use the libGraph graph library
#define GUROBI
#define GEDLIB

#include "SimplePatterns.h"
#include "Algorithms/GED/GEDEvaluation.h"
#include "Algorithms/GED/GEDLIBWrapper.h"
#include "GraphDataStructures/GraphBase.h"
#include "src/env/ged_env.hpp"

// main function to run the examples

int main() {
    const GraphStruct triangle = SimplePatterns::Circle(3);
    const GraphStruct square = SimplePatterns::Circle(5);

    GraphData<GraphStruct> graph_data;
    graph_data.add(triangle);
    graph_data.add(square);

    ged::GEDEnv<ged::LabelID, ged::LabelID, ged::LabelID> env;
    env.set_edit_costs(ged::Options::EditCosts::CONSTANT);
    AddGraphsToGEDEnvironment(env, graph_data);
    env.init();
    env.set_method(ged::Options::GEDMethod::F2);
    env.init_method();
    constexpr int source_id = 0;
    constexpr int target_id = 1;

    env.run_method(source_id,target_id);
    const GEDEvaluation result = ComputeGEDResult(env, graph_data, source_id, target_id);
    ged::NodeMap node_map = env.get_node_map(0,1);
    env.compute_induced_cost(0, 1, node_map);
    std::cout << "Approximated Distance: " << result.distance << std::endl;
    std::cout << "Time: " << result.time << " seconds" << std::endl;
    // print node mapping
    std::cout << "Node Mapping First: " << std::endl;
    for (const auto& x : result.node_mapping.first) {
        std::cout << x << " ";
    }
    std::cout << std::endl;
    std::cout << "Node Mapping Second: " << std::endl;
    for (const auto& x : result.node_mapping.second) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    EditPath edit_path;
    result.get_edit_path(edit_path, 0);
    // print edit path length
    std::cout << "Edit Path Length: " << edit_path.edit_path_graphs.size() - 1 << std::endl;
    GraphData<GraphStruct> edit_path_graphs;
    for (const auto& g : edit_path.edit_path_graphs) {
        edit_path_graphs.add(g);
    }
    edit_path_graphs.graphData.back().SetName(edit_path.target_graph.GetName());
    std::cout << edit_path << std::endl;
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

