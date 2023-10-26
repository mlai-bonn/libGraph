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
    GraphStruct graph = GraphStruct(10, {});

    graph.add_edge(0, 1);
    graph.add_edge(0, 2);
    graph.add_edge(8, 9);



    int counter = 0;
    std::vector<std::pair<NodeId, NodeId> > edges = {{0, 1}, {0, 2}, {8, 9}};
    GraphStruct::EdgeIterator edge = graph.first_edge();
    NodeId last_edge = graph.last_edge();
    EXPECT_EQ(edges[0], *edge);
    EXPECT_EQ(last_edge, 10);
    ++edge;
    EXPECT_EQ(edges[1], *edge);
    ++edge;
    EXPECT_EQ(edges[2], *edge);
    for (auto edge = graph.first_edge(); edge != graph.last_edge(); ++edge){
        EXPECT_EQ(edges[counter], *edge);
        ++counter;
    }
}

TEST(GraphIteratorsTestSuite, ExampleNodeIterator){
    GraphStruct graph = GraphStruct(10000, {});
    std::mt19937_64 generator(0);
    for (int i = 0; i < 10000; ++i) {
        for (int j = 0; j < 5000; ++j) {
            graph.add_edge(i, std::uniform_int_distribution<NodeId>(0, 9999)(generator));
        }
    }

    // iterator runtime for nodes
    auto start = std::chrono::high_resolution_clock::now();
    int counter = 0;
    for (auto node : graph){
        ++counter;
    }
    auto iterator_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "Iterator time for nodes: " << iterator_time << std::endl;
    EXPECT_EQ(counter, 10000);
    start = std::chrono::high_resolution_clock::now();
    counter = 0;
    for(int i = 0; i < graph.nodes(); ++i){
        ++counter;
    }
    auto for_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
    std::cout << "For time for nodes: " << for_time << std::endl;
    EXPECT_EQ(counter, 10000);
    std::cout << "Iterator time for nodes is " << for_time / iterator_time << " times faster than for loop" << std::endl;

}


#endif //GOOGLE_TESTS_GRAPHCLASSTEST_H
