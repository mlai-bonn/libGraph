//
// Created by florian on 24.10.23.
//

#ifndef GOOGLE_TESTS_GRAPHALGORITHMSTESTS_H
#define GOOGLE_TESTS_GRAPHALGORITHMSTESTS_H
#include <gtest/gtest.h>
#include "../include/libGraph.h"

TEST(SpanningSubgraphsTestSuite, ExampleBFSTreeUnconnected){
    GraphStruct graph = GraphStruct();
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    int components = BFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 0);
    EXPECT_EQ(tree.edges(), 0);
    EXPECT_NE(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }

    graph = GraphStruct(10, {});
    components = BFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 0);
    EXPECT_NE(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 10);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(distances[i], 0);
    }
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }

    graph = GraphStruct(4, {});
    graph.add_edge(0, 1);
    graph.add_edge(2, 3);
    components = BFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 4);
    EXPECT_EQ(tree.edges(), 2);
    EXPECT_NE(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 2);
    EXPECT_EQ(distances[0], 0);
    EXPECT_EQ(distances[1], 1);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }
}

TEST(SpanningSubgraphsTestSuite, ExampleRandomBFSTree){
    GraphStruct graph = SimplePatterns::ErdosRenyi(10, 20, 0, true);
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    int components = BFSSpanningTree(graph, tree, 0, visited, distances, false, 0);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 9);
    EXPECT_EQ(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 1);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }
}

TEST(SpanningSubgraphsTestSuite, ExampleDeterministicBFSTree){
    GraphStruct graph = SimplePatterns::ErdosRenyi(10, 20, 0, true);
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    int components = BFSSpanningTree(graph, tree, 0, visited, distances, true);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 9);
    EXPECT_EQ(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 1);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }
}

TEST(SpanningSubgraphsTestSuite, ExampleRandomDFSTree){
    GraphStruct graph = SimplePatterns::ErdosRenyi(10, 20, 0, true);
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    int components = DFSSpanningTree(graph, tree, 0, visited, distances, false, 0);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 9);
    EXPECT_EQ(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 1);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }
}

TEST(SpanningSubgraphsTestSuite, ExampleDeterministicDFSTree){
    GraphStruct graph = SimplePatterns::ErdosRenyi(10, 20, 0, true);
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    int components = DFSSpanningTree(graph, tree, 0, visited, distances, true);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 9);
    EXPECT_EQ(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 1);
    EXPECT_EQ(distances[0], 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }
}

TEST(SpanningSubgraphsTestSuite, ExampleDFSTreeUnconnected){
    GraphStruct graph = GraphStruct();
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    int components = BFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 0);
    EXPECT_EQ(tree.edges(), 0);
    EXPECT_NE(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }

    graph = GraphStruct(10, {});
    components = DFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 0);
    EXPECT_NE(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 10);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(distances[i], 0);
    }
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }

    graph = GraphStruct(4, {});
    graph.add_edge(0, 1);
    graph.add_edge(2, 3);
    components = DFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 4);
    EXPECT_EQ(tree.edges(), 2);
    EXPECT_NE(tree.get_type(), GraphType::TREE);
    EXPECT_EQ(components, 2);
    EXPECT_EQ(distances[0], 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }
}


TEST(BiconnectedComponentTestSuite, ExampleBiconnectedComponent){

    GraphStruct outer_planar_graph = Chepoi_Outerplanar_Graph();
    std::vector<std::vector<NodeId>> components;
    GetBiconnectedComponents(outer_planar_graph, components);
    EXPECT_EQ(components.size(), 1);
    EXPECT_EQ(components[0].size(), 14);
    GraphStruct bi_conn_wiki_graph = Bi_Conn_Wiki_Graph();
    GetBiconnectedComponents(bi_conn_wiki_graph, components);
    EXPECT_EQ(components.size(), 7);
    EXPECT_EQ(components[0].size(), 4);
    EXPECT_EQ(components[1].size(), 2);
    EXPECT_EQ(components[2].size(), 2);
    EXPECT_EQ(components[3].size(), 2);
    EXPECT_EQ(components[4].size(), 2);
    EXPECT_EQ(components[5].size(), 2);
    EXPECT_EQ(components[6].size(), 6);
}

TEST(OuterplanarSubgraphTestSuite, ExampleOuterplanarGraphData) {
    GraphStruct bi_conn_wiki_graph = Bi_Conn_Wiki_Graph();
    OuterplanarGraphData outerplanarGraphData = OuterplanarGraphData(bi_conn_wiki_graph);
    EXPECT_EQ(outerplanarGraphData.Components.size(), 2);
    EXPECT_EQ(outerplanarGraphData.get_bbTree().get_type(), GraphType::TREE);
}



#endif //GOOGLE_TESTS_GRAPHALGORITHMSTESTS_H
