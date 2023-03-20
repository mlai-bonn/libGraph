//
// Created by florian on 23.01.23.
//

#ifndef TESTGRAPHLIB_GRAPHFUNCTIONS_H
#define TESTGRAPHLIB_GRAPHFUNCTIONS_H

#include "../include/SimplePatterns.h"

bool ConversionTest(){
    for(std::string name : {"amazon", "dblp", "lj", "orkut", "youtube"}) {
        GraphStruct::Convert("../../../GraphData/Hops/com-" + name + ".ungraph.bin", "", GraphFormat::BGFS);
    }
    GraphStruct::Convert("../../../GraphData/Hops/", GraphFormat::BGFS);
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




#endif //TESTGRAPHLIB_GRAPHFUNCTIONS_H
