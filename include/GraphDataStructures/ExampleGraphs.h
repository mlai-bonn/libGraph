//
// Created by florian on 23.10.23.
//

#ifndef GOOGLE_TESTS_EXAMPLEGRAPHS_H
#define GOOGLE_TESTS_EXAMPLEGRAPHS_H



GraphStruct Bi_Conn_Wiki_Graph();
GraphStruct Chepoi_Outerplanar_Graph();

inline GraphStruct Bi_Conn_Wiki_Graph() {
    GraphStruct graph = GraphStruct(14, {});
    graph.add_edge(0, 1, false);

    graph.add_edge(0, 2, false);

    graph.add_edge(1, 3, false);

    graph.add_edge(3, 4, false);
    graph.add_edge(3, 5, false);
    graph.add_edge(4, 6, false);
    graph.add_edge(5, 6, false);

    graph.add_edge(2, 7, false);
    graph.add_edge(2, 8, false);
    graph.add_edge(2, 9, false);

    graph.add_edge(9, 10, false);
    graph.add_edge(9, 11, false);
    graph.add_edge(10, 11, false);
    graph.add_edge(10, 12, false);
    graph.add_edge(12, 8, false);
    graph.add_edge(8, 13, false);
    return graph;
}

inline GraphStruct Chepoi_Outerplanar_Graph() {
    GraphStruct outerplanar_graph = GraphStruct(15, {});
    outerplanar_graph.add_edge(0, 1, false);
    outerplanar_graph.add_edge(0, 2, false);
    outerplanar_graph.add_edge(1, 2, false);
    outerplanar_graph.add_edge(1, 3, false);
    outerplanar_graph.add_edge(2, 4, false);
    outerplanar_graph.add_edge(3, 5, false);
    outerplanar_graph.add_edge(4, 6, false);
    outerplanar_graph.add_edge(4, 7, false);
    outerplanar_graph.add_edge(5, 8, false);
    outerplanar_graph.add_edge(5, 9, false);
    outerplanar_graph.add_edge(6, 10, false);
    outerplanar_graph.add_edge(6, 11, false);
    outerplanar_graph.add_edge(7, 12, false);
    outerplanar_graph.add_edge(8, 13, false);
    outerplanar_graph.add_edge(9, 13, false);
    outerplanar_graph.add_edge(9, 10, false);
    outerplanar_graph.add_edge(11, 14, false);
    outerplanar_graph.add_edge(12, 14, false);
    return outerplanar_graph;
}


#endif //GOOGLE_TESTS_EXAMPLEGRAPHS_H
