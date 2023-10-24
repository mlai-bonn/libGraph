//
// Created by florian on 24.10.23.
//

#ifndef GOOGLE_TESTS_GRAPHCLOSURETESTS_H
#define GOOGLE_TESTS_GRAPHCLOSURETESTS_H

TEST(GraphClosureTestSuite, ExampleGraphClosureUnconnected){
    GraphStruct graph = GraphStruct(2,{});
    GraphClosureSP graphClosure = GraphClosureSP(graph);
    ClosureParameters closureParameters = {{0}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 1);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0}));

    closureParameters = {{0, 1}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 2);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, 1}));
}

TEST(GraphClosureTestSuite, ExampleGraphClosureConnectedSimple){
    GraphStruct graph = GraphStruct(2,{});
    graph.add_edge(0, 1);
    GraphClosureSP graphClosure = GraphClosureSP(graph);
    ClosureParameters closureParameters = {{0}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 1);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0}));

    closureParameters = {{0, 1}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 2);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, 1}));
}

TEST(GraphClosureTestSuite, ExampleGraphClosureConnectedAdvanced){
    GraphStruct graph = SimplePatterns::Path(50);
    GraphClosureSP graphClosure = GraphClosureSP(graph);
    ClosureParameters closureParameters = {{0}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 1);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0}));

    closureParameters = {{0, 49}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 50);
    closureParameters = {{1, 3}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 3);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({1, 2, 3}));

    graph = SimplePatterns::Circle(4);
    closureParameters = {{0, 2}};
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 4);


    for (int i = 3; i < 20; ++i) {
        graph = SimplePatterns::Circle(i);
        closureParameters = {{0, (NodeId) i/2 }};
        graphClosure.closure(closureParameters);
        // if i is odd
        if (i % 2 == 1) {
            EXPECT_EQ(closureParameters.closed_set.size(), (size_t) i/2 + 1);
        }
        else{
            EXPECT_EQ(closureParameters.closed_set.size(), i);
        }

        graph = SimplePatterns::StarGraph(i);
        closureParameters = {{0, (NodeId) i-1}};
        graphClosure.closure(closureParameters);
        EXPECT_EQ(closureParameters.closed_set.size(), 2);
        EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, (NodeId) i-1}));
        closureParameters = {{1, (NodeId) i-1}};
        graphClosure.closure(closureParameters);
        EXPECT_EQ(closureParameters.closed_set.size(), 3);
        EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, 1, (NodeId) i-1}));

        graph = SimplePatterns::MaxPreClosure(i);
        for (int j = 0; j < graph.nodes(); ++j) {
            for (int k = 0; k < graph.nodes(); ++k) {
                closureParameters = {{(NodeId) j, (NodeId) k}};
                graphClosure.closure(closureParameters);
                if (j==k){
                    EXPECT_EQ(closureParameters.closed_set.size(), 1);
                    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({(NodeId) j}));
                }
                else if (graph.edge(j, k)){
                    EXPECT_EQ(closureParameters.closed_set.size(), 2);
                    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({(NodeId) j, (NodeId) k}));
                }
                else{
                    EXPECT_EQ(closureParameters.closed_set.size(), graph.nodes());
                }
            }
        }
    }




}

#endif //GOOGLE_TESTS_GRAPHCLOSURETESTS_H
