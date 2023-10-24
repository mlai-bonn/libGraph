//
// Created by florian on 23.10.23.
//

#ifndef GOOGLE_TESTS_EXAMPLEGRAPHS_H
#define GOOGLE_TESTS_EXAMPLEGRAPHS_H



GraphStruct Bi_Conn_Wiki_Graph();
GraphStruct Chepoi_Outerplanar_Graph();

inline GraphStruct Bi_Conn_Wiki_Graph() {
    GraphStruct graph = GraphStruct(14, {});
    graph.add_edge_no_check(0, 1);

    graph.add_edge_no_check(0, 2);

    graph.add_edge_no_check(1, 3);

    graph.add_edge_no_check(3, 4);
    graph.add_edge_no_check(3, 5);
    graph.add_edge_no_check(4, 6);
    graph.add_edge_no_check(5, 6);

    graph.add_edge_no_check(2, 7);
    graph.add_edge_no_check(2, 8);
    graph.add_edge_no_check(2, 9);

    graph.add_edge_no_check(9, 10);
    graph.add_edge_no_check(9, 11);
    graph.add_edge_no_check(10, 11);
    graph.add_edge_no_check(10, 12);
    graph.add_edge_no_check(12, 8);
    graph.add_edge_no_check(8, 13);
    return graph;
}

inline GraphStruct Chepoi_Outerplanar_Graph() {
    GraphStruct outerplanar_graph = GraphStruct(15, {});
    outerplanar_graph.add_edge_no_check(0, 1);
    outerplanar_graph.add_edge_no_check(0, 2);
    outerplanar_graph.add_edge_no_check(1, 2);
    outerplanar_graph.add_edge_no_check(1, 3);
    outerplanar_graph.add_edge_no_check(2, 4);
    outerplanar_graph.add_edge_no_check(3, 5);
    outerplanar_graph.add_edge_no_check(4, 6);
    outerplanar_graph.add_edge_no_check(4, 7);
    outerplanar_graph.add_edge_no_check(5, 8);
    outerplanar_graph.add_edge_no_check(5, 9);
    outerplanar_graph.add_edge_no_check(6, 10);
    outerplanar_graph.add_edge_no_check(6, 11);
    outerplanar_graph.add_edge_no_check(7, 12);
    outerplanar_graph.add_edge_no_check(8, 13);
    outerplanar_graph.add_edge_no_check(9, 13);
    outerplanar_graph.add_edge_no_check(9, 10);
    outerplanar_graph.add_edge_no_check(11, 14);
    outerplanar_graph.add_edge_no_check(12, 14);
    return outerplanar_graph;
}


#endif //GOOGLE_TESTS_EXAMPLEGRAPHS_H
