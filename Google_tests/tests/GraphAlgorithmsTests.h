//
// Created by florian on 24.10.23.
//

#ifndef GOOGLE_TESTS_GRAPHALGORITHMSTESTS_H
#define GOOGLE_TESTS_GRAPHALGORITHMSTESTS_H
#include <gtest/gtest.h>
#include "libGraph.h"

TEST(SpanningSubgraphsTestSuite, ExampleBFSTreeUnconnected){
    GraphStruct graph = GraphStruct();
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    int components = BFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 0);
    EXPECT_EQ(tree.edges(), 0);
    EXPECT_NE(tree.GetType(), GraphType::TREE);
    EXPECT_EQ(components, 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }

    graph = GraphStruct(10, {});
    components = BFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 0);
    EXPECT_NE(tree.GetType(), GraphType::TREE);
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
    EXPECT_NE(tree.GetType(), GraphType::TREE);
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
    EXPECT_EQ(tree.GetType(), GraphType::TREE);
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
    EXPECT_EQ(tree.GetType(), GraphType::TREE);
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
    EXPECT_EQ(tree.GetType(), GraphType::TREE);
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
    EXPECT_EQ(tree.GetType(), GraphType::TREE);
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
    EXPECT_NE(tree.GetType(), GraphType::TREE);
    EXPECT_EQ(components, 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }

    graph = GraphStruct(10, {});
    components = DFSSpanningTree(graph, tree, 0, visited, distances);
    EXPECT_EQ(tree.nodes(), 10);
    EXPECT_EQ(tree.edges(), 0);
    EXPECT_NE(tree.GetType(), GraphType::TREE);
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
    EXPECT_NE(tree.GetType(), GraphType::TREE);
    EXPECT_EQ(components, 2);
    EXPECT_EQ(distances[0], 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        EXPECT_EQ(visited[i], true);
    }
}


TEST(BiconnectedComponentTestSuite, ExampleBiconnectedComponent){

    GraphStruct outer_planar_graph = Chepoi_Outerplanar_Graph();
    std::vector<std::vector<NodeId>> components;
    GraphAlgorithms::GetBiconnectedComponents(outer_planar_graph, components);
    EXPECT_EQ(components.size(), 1);
    EXPECT_EQ(components[0].size(), 15);
    GraphStruct bi_conn_wiki_graph = Bi_Conn_Wiki_Graph();
    GraphAlgorithms::GetBiconnectedComponents(bi_conn_wiki_graph, components);
    EXPECT_EQ(components.size(), 7);
    EXPECT_EQ(components[0].size(), 4);
    EXPECT_EQ(components[1].size(), 2);
    EXPECT_EQ(components[2].size(), 2);
    EXPECT_EQ(components[3].size(), 2);
    EXPECT_EQ(components[4].size(), 6);
    EXPECT_EQ(components[5].size(), 2);
    EXPECT_EQ(components[6].size(), 2);
}

TEST(BiconnectedComponentTestSuite, ExampleChainOfCircles) {
    int size = 3;
    int circle_number = 3;
    GraphStruct n_circles = GraphStruct(circle_number * (size - 1) + 1, {});
    for (int i = 0; i < circle_number; ++i) {
        for (int j = 0; j < size; ++j) {
            n_circles.add_edge(i * (size - 1) + j, i * (size - 1) + (j + 1) % size);
        }
    }
    std::vector<std::vector<NodeId>> components;
    GraphAlgorithms::GetBiconnectedComponents(n_circles, components);
    EXPECT_EQ(components.size(), circle_number);
    for (int i = 0; i < circle_number; ++i) {
        EXPECT_EQ(components[i].size(), size);
    }
    std::vector<GraphStruct> graph_components;
    GraphAlgorithms::GetBiconnectedComponents(n_circles, graph_components);
    EXPECT_EQ(graph_components.size(), circle_number);
    for (int i = 0; i < circle_number; ++i) {
        EXPECT_EQ(graph_components[i].nodes(), size);
        EXPECT_EQ(graph_components[i].edges(), size);
    }
}

TEST(BiconnectedComponentTestSuite, ExampleCircleWithSpikes) {
    int size = 3;
    int spikes = 3;
    GraphStruct spike_circle = SimplePatterns::Circle(size);
    for (int i = 0; i < std::min(size, spikes); ++i) {
        NodeId new_node = spike_circle.add_node();
        spike_circle.add_edge(i, new_node);
    }
    std::vector<std::vector<NodeId>> components;
    GraphAlgorithms::GetBiconnectedComponents(spike_circle, components);
    EXPECT_EQ(components.size(), size + 1);
    EXPECT_EQ(components[0].size(), 2);
    EXPECT_EQ(components[1].size(), 2);
    EXPECT_EQ(components[2].size(), 2);
    EXPECT_EQ(components[3].size(), 3);
    std::vector<GraphStruct> graph_components;
    GraphAlgorithms::GetBiconnectedComponents(spike_circle, graph_components);
    EXPECT_EQ(graph_components.size(), size + 1);
    EXPECT_EQ(graph_components[0].nodes(), 2);
    EXPECT_EQ(graph_components[0].edges(), 1);
    EXPECT_EQ(graph_components[1].nodes(), 2);
    EXPECT_EQ(graph_components[1].edges(), 1);
    EXPECT_EQ(graph_components[2].nodes(), 2);
    EXPECT_EQ(graph_components[2].edges(), 1);
    EXPECT_EQ(graph_components[3].nodes(), 3);
    EXPECT_EQ(graph_components[3].edges(), 3);
}

TEST(BiconnectedComponentTestSuite, ExampleCircleWithSpikesGraphExtended) {
    int size = 3;
    int spikes = 3;
    GraphExtended spike_circle = GraphExtended(SimplePatterns::Circle(size));
    for (int i = 0; i < std::min(size, spikes); ++i) {
        NodeId new_node = spike_circle.add_node();
        spike_circle.add_edge(i, new_node);
    }
    std::vector<std::vector<NodeId>> components;
    GraphAlgorithms::GetBiconnectedComponents(spike_circle, components);
    EXPECT_EQ(components.size(), size + 1);
    EXPECT_EQ(components[0].size(), 2);
    EXPECT_EQ(components[1].size(), 2);
    EXPECT_EQ(components[2].size(), 2);
    EXPECT_EQ(components[3].size(), 3);
    std::vector<GraphStruct> graph_components;
    GraphAlgorithms::GetBiconnectedComponents(spike_circle, graph_components);
    EXPECT_EQ(graph_components.size(), size + 1);
    EXPECT_EQ(graph_components[0].nodes(), 2);
    EXPECT_EQ(graph_components[0].edges(), 1);
    EXPECT_EQ(graph_components[1].nodes(), 2);
    EXPECT_EQ(graph_components[1].edges(), 1);
    EXPECT_EQ(graph_components[2].nodes(), 2);
    EXPECT_EQ(graph_components[2].edges(), 1);
    EXPECT_EQ(graph_components[3].nodes(), 3);
    EXPECT_EQ(graph_components[3].edges(), 3);
}

TEST(BiconnectedComponentTestSuite, ExampleCirclesWithPathConnections) {
    int size = 3;
    int number = 2;
    int path_length = 3;
    GraphStruct graphStruct = GraphStruct(size*number + (path_length-1)*(number - 1), {});

    // add circle edges
    for (int i = 0; i < number; ++i) {
        for (int j = 0; j < size; ++j) {
            NodeId src = i * size + i * (path_length - 1) + j;
            NodeId dst = i * size + i * (path_length - 1) + (j + 1) % size;
            graphStruct.add_edge(src, dst);
        }
    }
    // add path edges
    for (int i = 0; i < number - 1; ++i) {
        for (int j = 0; j < path_length; ++j) {
            NodeId src = i * (size + path_length - 1) + size-1 + j;
            NodeId dst = src + 1;
            graphStruct.add_edge(src, dst);
        }
    }
    std::vector<std::vector<NodeId>> components;
    GraphAlgorithms::GetBiconnectedComponents(graphStruct, components);
    EXPECT_EQ(components.size(), number + (number - 1) * path_length);
    std::vector<GraphStruct> graph_components;
    GraphAlgorithms::GetBiconnectedComponents(graphStruct, graph_components);
    EXPECT_EQ(graph_components.size(), number + (number - 1) * path_length);
}

TEST(BiconnectedComponentTestSuite, ExampleCirclesWithPathConnectionsGraphExtended) {
    int size = 3;
    int number = 2;
    int path_length = 3;
    GraphExtended graphStruct = GraphExtended(size*number + (path_length-1)*(number - 1), {});

    // add circle edges
    for (int i = 0; i < number; ++i) {
        for (int j = 0; j < size; ++j) {
            NodeId src = i * size + i * (path_length - 1) + j;
            NodeId dst = i * size + i * (path_length - 1) + (j + 1) % size;
            graphStruct.add_edge(src, dst);
        }
    }
    // add path edges
    for (int i = 0; i < number - 1; ++i) {
        for (int j = 0; j < path_length; ++j) {
            NodeId src = i * (size + path_length - 1) + size-1 + j;
            NodeId dst = src + 1;
            graphStruct.add_edge(src, dst);
        }
    }
    std::vector<std::vector<NodeId>> components;
    GraphAlgorithms::GetBiconnectedComponents(graphStruct, components);
    EXPECT_EQ(components.size(), number + (number - 1) * path_length);
    std::vector<GraphStruct> graph_components;
    GraphAlgorithms::GetBiconnectedComponents(graphStruct, graph_components);
    EXPECT_EQ(graph_components.size(), number + (number - 1) * path_length);
}

TEST(OuterplanarSubgraphTestSuite, ExampleOuterplanarGraphData) {
    GraphStruct bi_conn_wiki_graph = Bi_Conn_Wiki_Graph();
    OuterplanarGraphData outerplanarGraphData = OuterplanarGraphData(bi_conn_wiki_graph);
    EXPECT_EQ(outerplanarGraphData.Components.size(), 2);
    EXPECT_EQ(outerplanarGraphData.get_bbTree().GetType(), GraphType::TREE);
}



#endif //GOOGLE_TESTS_GRAPHALGORITHMSTESTS_H
