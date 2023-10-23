//
// Created by florian on 23.01.23.
//

#ifndef TESTGRAPHLIB_GRAPHFUNCTIONS_H
#define TESTGRAPHLIB_GRAPHFUNCTIONS_H

#include "../../include/SimplePatterns.h"
#include "../../include/io/io.h"
#include "../../include/io/FileEvaluation.h"
#include "../../include/Algorithms/GraphMetricApproximation/EvaluationGraphMetricApproximation.h"

bool ConversionTest(){
    for(std::string name : {"amazon", "dblp", "lj", "orkut", "youtube"}) {
        GraphStruct::Convert("../../../GraphData/Hops/com-" + name + ".ungraph.bin", "", GraphFormat::BGFS);
    }
    GraphStruct::Convert("../../../GraphData/Hops/", GraphFormat::BGFS);
    return true;
}

bool BGFS_to_TXT(){
    for(auto i : std::vector<int>({4, 5, 6, 7, 8, 9}))
    {
        std::string name = "size" + std::to_string(i);
        GraphData graphs = GraphData<GraphStruct>("../../../../GraphData/Hops/patterns/" + name + ".bgfs");
        int counter = 0;
        std::filesystem::create_directory("../../../../hops/Patterns" + std::to_string(i) + "/");
        for(auto& g : graphs.graphData){
            if (g.IsTree()) {
                SaveParams params = {"../../../../hops/Patterns" + std::to_string(i) + "/", "tree_" + std::to_string(counter), GraphFormat::EDGES};
                g.Save(params);
                ++counter;
            }
        }
    }

    return true;
}

bool LargestComponentTest(){
    GraphStruct graphStruct = GraphStruct("../../../GraphData/Hops/com-orkut.ungraph.txt");
    auto largestComponent = GraphStruct::GetLargestComponent(graphStruct);
    largestComponent.Save({});
    return true;
}

double edge_weight(const DDataGraph& graph, EDGE edge){
    return graph.get_edge_data(edge, "weight");
}

bool DijkstraTest() {
    DDataGraph dGraph = DDataGraph("../test/TestData/DGraph.edges",true,"../test/TestData/DGraph.data");
    std::vector<double> distances = std::vector<double>(dGraph.nodes(), 0);
    bool connected = DDataGraph::dijkstra(dGraph, 3, &edge_weight, distances);
    std::vector<NodeId> path;
    double length;
    bool found_path = DDataGraph::dijkstra(dGraph, 2, 4, &edge_weight, path, distances, length);
    return connected && found_path;
}

bool ErdosRenyi(){
    int nodes = 20, edges = 100;
    int num = 10;
    GraphData<GraphStruct> graphData;
    for (int i = 0; i < num; ++i) {
        GraphStruct graphStruct = SimplePatterns::ErdosRenyi(nodes,edges,i);
        graphData.add(graphStruct);
    }
    graphData.Save({"../../../GraphData/ErdosRenyi/","er_" + std::to_string(nodes) + "_" + std::to_string(edges), GraphFormat::BGFS});
    return true;
}

bool DFS(){
    GraphStruct graph = GraphStruct("../../../../GraphData/Hops/com-amazon.ungraph.bgfs");
    GraphStruct tree;
    Nodes order;
    GraphStruct::DFS(graph, tree, order);
    int x = 0;
}

bool MapGraph(){

    GraphStruct g(5,{});
    g.add_edge(0,1);
    g.add_edge(2,0);
    g.add_edge(0, 3);
    g.add_edge(3,2);
    g.add_edge(4, 2);

    GraphStruct t;
    Nodes o;
    GraphStruct::DFS(g, t, o);
    GraphStruct::ReorderGraph(g, {4,2,0,1,3});


    GraphStruct graph = GraphStruct("../../../../GraphData/Hops/com-amazon.ungraph.bgfs");
    Nodes order;
    GraphStruct tree;
    GraphStruct::DFS(graph, tree, order);
    Nodes newOrder = Nodes(order.size());
    int i = 0;
    for (auto o : order) {
        newOrder[o] = i;
        ++i;
    }
    GraphStruct newGraph(graph);
    GraphStruct::ReorderGraph(newGraph, order);
    GraphStruct::ReorderGraph(newGraph, newOrder);

    bool y = newGraph == graph;

    int x = 0;
}

bool GenerateGraphs(){
    for (int size : {100, 1000, 10000}) {
        GraphStruct g = GraphStruct(SimplePatterns::FullyConnected(size));
        g.Save({"../../../../GraphData/Hops/", "fc_" + std::to_string(size), GraphFormat::BGFS});
  }
    return true;
}

bool GraphsToLatex(){

    GraphData graphs = GraphData<GraphStruct>(GraphFormat::BGFS, "../../../../GraphData/Hops/", "");

    FileEvaluation evaluation = FileEvaluation("../../../../ChoPS/final_results/", "graphs");
    //std::vector<std::vector<std::string>> info;

    //info.emplace_back(std::vector<std::string>{"Name", "Size", "Edges", "Avg. Degree", "Max. Degree"});

    for (auto& graph : graphs.graphData) {

        evaluation.headerValueInsert({"Name", "Size", "Edges", "Avg. Degree", "Max. Degree"},
                                     {graph.GetName(), std::to_string(graph.nodes()), std::to_string(graph.edges()),
                                      std::to_string((graph.edges())*2/graph.nodes()),
                                      std::to_string(graph.maxDegree)});

//        std::vector<std::string> g_info;
//        g_info.emplace_back(graph.GetName());
//        g_info.emplace_back(std::to_string(graph.nodes()));
//        g_info.emplace_back(std::to_string(graph.edges()));
//        g_info.emplace_back(std::to_string((graph.edges())*2/graph.nodes()));
//        g_info.emplace_back(std::to_string(graph.maxDegree));
//        info.emplace_back(g_info);
    }
    evaluation.save();
    evaluation.to_latex_table();
    //ToLatexTable("../../../../ChoPS/final_results/graphs.tex", info);
}

bool TestLayeringTree(){
    // generate some random graph
    // GraphStruct graph = SimplePatterns::ErdosRenyi(10, 20, 0, true);
    GraphStruct graph = GraphStruct(15, {});
    graph.add_edge(0,1);
    graph.add_edge(0,2);
    graph.add_edge(1,2);
    graph.add_edge(1,3);
    graph.add_edge(2,4);
    graph.add_edge(3,5);
    graph.add_edge(4,6);
    graph.add_edge(4,7);
    graph.add_edge(5,8);
    graph.add_edge(5,9);
    graph.add_edge(6,10);
    graph.add_edge(6,11);
    graph.add_edge(7,12);
    graph.add_edge(8,13);
    graph.add_edge(9,13);
    graph.add_edge(9,10);
    graph.add_edge(11,14);
    graph.add_edge(12,14);

    // Fully connected graph
    //graph = SimplePatterns::ErdosRenyi(10000, 100000, 0, true);

    // print graph has been generated
    std::cout << "Graph has been generated" << std::endl;

    EvaluateApproximations evaluateApproximations(graph);

    for (int i = 0; i < graph.nodes(); ++i) {
        ParametersEvaluation params({ApproximationType::LAYERING_TREE, 0, i, 1, false});
        evaluateApproximations.Evaluation(params);
    }

    ParametersEvaluation parametersEvaluation({ApproximationType::RANDOM_SPANNING_TREES, 0, 0, 100, false});
    evaluateApproximations.Evaluation(parametersEvaluation);
    return true;
}

bool TestSpanningTrees(){
    // generate some random graph
    // GraphStruct graph = SimplePatterns::ErdosRenyi(10, 20, 0, true);
    GraphStruct graph = GraphStruct(15, {});
    graph.add_edge(0,1);
    graph.add_edge(0,2);
    graph.add_edge(1,2);
    graph.add_edge(1,3);
    graph.add_edge(2,4);
    graph.add_edge(3,5);
    graph.add_edge(4,6);
    graph.add_edge(4,7);
    graph.add_edge(5,8);
    graph.add_edge(5,9);
    graph.add_edge(6,10);
    graph.add_edge(6,11);
    graph.add_edge(7,12);
    graph.add_edge(8,13);
    graph.add_edge(9,13);
    graph.add_edge(9,10);
    graph.add_edge(11,14);
    graph.add_edge(12,14);

    // TODO
    return true;
}

bool TestGraphAlgorithms(){
    // generate some random graph
    // GraphStruct graph = SimplePatterns::ErdosRenyi(10, 20, 0, true);
    GraphStruct graph = GraphStruct(15, {});
    graph.add_edge(0,1);
    graph.add_edge(0,2);
    graph.add_edge(1,2);
    graph.add_edge(1,3);
    graph.add_edge(2,4);
    graph.add_edge(3,5);
    graph.add_edge(4,6);
    graph.add_edge(4,7);
    graph.add_edge(5,8);
    graph.add_edge(5,9);
    graph.add_edge(6,10);
    graph.add_edge(6,11);
    graph.add_edge(7,12);
    graph.add_edge(8,13);
    graph.add_edge(9,13);
    graph.add_edge(9,10);
    graph.add_edge(11,14);
    graph.add_edge(12,14);



    // TODO
    return true;
}





#endif //TESTGRAPHLIB_GRAPHFUNCTIONS_H
