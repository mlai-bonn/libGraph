//
// Created by florian on 23.01.23.
//

#ifndef TESTGRAPHLIB_GRAPHFUNCTIONS_H
#define TESTGRAPHLIB_GRAPHFUNCTIONS_H

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


#endif //TESTGRAPHLIB_GRAPHFUNCTIONS_H
