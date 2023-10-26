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


#endif //GOOGLE_TESTS_GRAPHCLASSTEST_H
