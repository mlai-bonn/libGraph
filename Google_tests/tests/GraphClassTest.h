//
// Created by florian on 24.10.23.
//

#ifndef GOOGLE_TESTS_GRAPHCLASSTEST_H
#define GOOGLE_TESTS_GRAPHCLASSTEST_H

#include <gtest/gtest.h>
#include "GraphDataStructures/GraphBase.h"

TEST(GraphConstructorTestSuite, ExampleDefaultGraphConstructor){
    GraphStruct graph = GraphStruct();

    EXPECT_EQ(graph.nodes(), 0);
    EXPECT_EQ(graph.edges(), 0);
}

TEST(GraphConstructorTestSuite, ExampleGraphConstructor){
    GraphStruct graph = GraphStruct(10, {});

    EXPECT_EQ(graph.nodes(), 10);
    EXPECT_EQ(graph.edges(), 0);
}

TEST(GraphNodesTestSuite, ExampleAddNode){
    GraphStruct graph = GraphStruct(10, {});

    graph.add_node();

    EXPECT_EQ(graph.nodes(), 11);
    EXPECT_EQ(graph.edges(), 0);
}

TEST(GraphNodesTestSuite, ExampleAddNodes){
    GraphStruct graph = GraphStruct(10, {});

    graph.add_node(10);

    EXPECT_EQ(graph.nodes(), 20);
    EXPECT_EQ(graph.edges(), 0);
}

TEST(GraphEdgesTestSuite, ExampleAddEdge){
    GraphStruct graph = GraphStruct(10, {});

    graph.add_edge(0, 1);

    EXPECT_EQ(graph.nodes(), 10);
    EXPECT_EQ(graph.edges(), 1);
    EXPECT_EQ(graph.edge(0, 1), true);
    EXPECT_EQ(graph.degree(0), 1);
    EXPECT_EQ(graph.degree(1), 1);
    for (int i = 2; i < graph.nodes(); ++i) {
        EXPECT_EQ(graph.degree(i), 0);
    }
}
TEST(GraphEdgesTestSuite, ExampleAddEdgeNoCheck){
    GraphStruct graph = GraphStruct(10, {});

    graph.add_edge_no_check(0, 1);

    EXPECT_EQ(graph.nodes(), 10);
    EXPECT_EQ(graph.edges(), 1);
    EXPECT_EQ(graph.edge(0, 1), true);
    EXPECT_EQ(graph.degree(0), 1);
    EXPECT_EQ(graph.degree(1), 1);
    for (int i = 2; i < graph.nodes(); ++i) {
        EXPECT_EQ(graph.degree(i), 0);
    }
}

TEST(GraphEdgesTestSuite, ExampleRemoveEdge){
    GraphStruct graph = GraphStruct(10, {});

    graph.add_edge(0, 1);
    graph.remove_edge(0, 1);

    EXPECT_EQ(graph.nodes(), 10);
    EXPECT_EQ(graph.edges(), 0);
    EXPECT_EQ(graph.edge(0, 1), false);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(graph.degree(i), 0);
    }
}

TEST(GraphIteratorsTestSuite, ExampleIterateOverNodes){
    GraphStruct graph = GraphStruct(10, {});

    graph.add_edge(0, 1);

    int counter = 0;
    for (auto node : graph){
        EXPECT_EQ(node, counter);
        ++counter;
    }
}

TEST(GraphIteratorsTestSuite, ExampleEdgeIterator){
    GraphStruct no_edges = GraphStruct(100, {});
    int counter = 0;
    // Iterate over edges
    for (auto edge = no_edges.first_edge(); edge != no_edges.last_edge(); ++edge) {
        ++counter;
    }

    GraphStruct graph = GraphStruct(10, {});

    graph.add_edge(0, 1);
    graph.add_edge(0, 2);
    graph.add_edge(8, 9);



    counter = 0;
    std::vector<std::pair<NodeId, NodeId> > edges = {{0, 1}, {0, 2}, {8, 9}};
    GraphStruct::EdgeIterator edge = graph.first_edge();
    NodeId last_edge = graph.last_edge();
    EXPECT_EQ(edges[0], *edge);
    EXPECT_EQ(last_edge, 10);
    ++edge;
    EXPECT_EQ(edges[1], *edge);
    ++edge;
    EXPECT_EQ(edges[2], *edge);
    ++edge;
    for (auto edge = graph.first_edge(); edge != graph.last_edge(); ++edge){
        EXPECT_EQ(edges[counter], *edge);
        ++counter;
    }
}

TEST(GraphIteratorsTestSuite, ExampleRuntimeIterator){
    INDEX graph_size = 100000000;
    GraphStruct graph = GraphStruct(graph_size, {});

    // iterator runtime for nodes
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<NodeId> nodes;
    for (auto node : graph){
        nodes.emplace_back(node);
    }
    auto iterator_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Iterator time for nodes: " << ((double) iterator_time / 1e9) << std::endl;
    EXPECT_EQ(nodes.size(), graph_size);
    start = std::chrono::high_resolution_clock::now();
    nodes.clear();
    for(NodeId i = 0; i < graph.nodes(); ++i){
        nodes.emplace_back(i);
    }
    auto for_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "For time for nodes: " << ((double) for_time / 1e9) << std::endl;
    EXPECT_EQ(nodes.size(), graph_size);
    std::cout << "For loop is " << (double) iterator_time / (double) for_time << " times faster than the iterator" << std::endl;

    std::mt19937_64 generator(0);
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 100; ++j) {
            graph.add_edge(i, std::uniform_int_distribution<NodeId>(0, graph_size - 1)(generator));
        }
    }
    EXPECT_EQ(graph.edges(), 100000);

    std::vector<std::pair<NodeId, NodeId>> edges;
    start = std::chrono::high_resolution_clock::now();
    for (NodeId i = 0; i < graph.nodes(); ++i) {
        for (int j = 0; j < graph.degree(i); ++j) {
            NodeId neighbor = graph[i][j];
            edges.emplace_back(i, neighbor);
        }
    }
    auto edge_for_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Num Edges:" << edges.size() << std::endl;
    EXPECT_EQ(graph.edges()*2, edges.size());

    edges.clear();
    start = std::chrono::high_resolution_clock::now();
    for (auto i : graph) {
        for (int j = 0; j < graph.degree(i); ++j) {
            NodeId neighbor = graph[i][j];
            edges.emplace_back(i, neighbor);
        }
    }
    auto edge_node_iterator = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Num Edges:" << edges.size() << std::endl;
    EXPECT_EQ(graph.edges()*2, edges.size());

    edges.clear();
    start = std::chrono::high_resolution_clock::now();
    for (NodeId i = 0; i < graph.nodes(); ++i) {
        for (auto neighbor = graph.beginN(i); neighbor != graph.endN(i);++neighbor) {
            NodeId n = *neighbor;
            edges.emplace_back(i, n);
        }
    }
    auto neighbor_iterator_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Num Edges:" << edges.size() << std::endl;
    EXPECT_EQ(graph.edges()*2, edges.size());

    edges.clear();
    start = std::chrono::high_resolution_clock::now();
    for (auto edge = graph.first_edge(); edge != graph.last_edge(); ++edge) {
        edges.emplace_back(edge.srcId, edge.dstId);
    }
    auto edge_iterator_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Num Edges:" << edges.size() << std::endl;
    EXPECT_EQ(graph.edges(), edges.size());

    std::cout << "Iterating over all edges (for loop) is: " << (double) edge_for_time / (double) 1e9 << std::endl;
    std::cout << "For time for iterating over all edges is: " << (double) edge_for_time / (double) 1e9 << std::endl;
    std::cout << "For time for iterating over all edges using the neighbor iterator is: " << (double) neighbor_iterator_time / (double) 1e9 << std::endl;
    std::cout << "For time for iterating over all edges using the edge iterator is: " << (double) edge_iterator_time / (double) 1e9 << std::endl;

    // NeighborIterators


}


#endif //GOOGLE_TESTS_GRAPHCLASSTEST_H
