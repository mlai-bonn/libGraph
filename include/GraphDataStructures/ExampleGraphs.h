//
// Created by florian on 23.10.23.
//

#ifndef GOOGLE_TESTS_EXAMPLEGRAPHS_H
#define GOOGLE_TESTS_EXAMPLEGRAPHS_H



GraphStruct Bi_Conn_Wiki_Graph();
GraphStruct Chepoi_Outerplanar_Graph();

inline GraphStruct Bi_Conn_Wiki_Graph() {
    GraphStruct graph = GraphStruct(14, {});
    graph.AddEdge(0, 1, false);

    graph.AddEdge(0, 2, false);

    graph.AddEdge(1, 3, false);

    graph.AddEdge(3, 4, false);
    graph.AddEdge(3, 5, false);
    graph.AddEdge(4, 6, false);
    graph.AddEdge(5, 6, false);

    graph.AddEdge(2, 7, false);
    graph.AddEdge(2, 8, false);
    graph.AddEdge(2, 9, false);

    graph.AddEdge(9, 10, false);
    graph.AddEdge(9, 11, false);
    graph.AddEdge(10, 11, false);
    graph.AddEdge(10, 12, false);
    graph.AddEdge(12, 8, false);
    graph.AddEdge(8, 13, false);
    return graph;
}

inline GraphStruct Chepoi_Outerplanar_Graph() {
    GraphStruct outerplanar_graph = GraphStruct(15, {});
    outerplanar_graph.AddEdge(0, 1, false);
    outerplanar_graph.AddEdge(0, 2, false);
    outerplanar_graph.AddEdge(1, 2, false);
    outerplanar_graph.AddEdge(1, 3, false);
    outerplanar_graph.AddEdge(2, 4, false);
    outerplanar_graph.AddEdge(3, 5, false);
    outerplanar_graph.AddEdge(4, 6, false);
    outerplanar_graph.AddEdge(4, 7, false);
    outerplanar_graph.AddEdge(5, 8, false);
    outerplanar_graph.AddEdge(5, 9, false);
    outerplanar_graph.AddEdge(6, 10, false);
    outerplanar_graph.AddEdge(6, 11, false);
    outerplanar_graph.AddEdge(7, 12, false);
    outerplanar_graph.AddEdge(8, 13, false);
    outerplanar_graph.AddEdge(9, 13, false);
    outerplanar_graph.AddEdge(9, 10, false);
    outerplanar_graph.AddEdge(11, 14, false);
    outerplanar_graph.AddEdge(12, 14, false);
    return outerplanar_graph;
}


#endif //GOOGLE_TESTS_EXAMPLEGRAPHS_H
