//
// Created by florian on 24.10.23.
//

#ifndef GOOGLE_TESTS_GRAPHCLOSURETESTS_H
#define GOOGLE_TESTS_GRAPHCLOSURETESTS_H

TEST(GraphClosureTestSuite, ExampleGraphClosureUnconnected){
    GraphStruct graph = GraphStruct(2,{});
    GraphClosure graphClosure = GraphClosure(graph);
    GraphClosureParameters closureParameters = GraphClosureParameters({.input_set = {0}});
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 1);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0}));

    closureParameters = GraphClosureParameters({.input_set = {0, 1}});
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 2);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, 1}));
}

TEST(GraphClosureTestSuite, ExampleGraphClosureConnectedSimple){
    GraphStruct graph = GraphStruct(2,{});
    graph.add_edge(0, 1);
    GraphClosure graphClosure = GraphClosure(graph);
    GraphClosureParameters closureParameters = GraphClosureParameters({.input_set = {0}});
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 1);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0}));

    closureParameters =  GraphClosureParameters({.input_set ={0, 1}});
    graphClosure.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 2);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, 1}));
}

TEST(GraphClosureTestSuite, ExampleGraphClosureConnectedAdvanced){
    GraphStruct path = SimplePatterns::Path(50);
    GraphClosure graphClosurePath = GraphClosure(path);
    GraphClosureParameters closureParameters = GraphClosureParameters({.input_set = {0}});
    graphClosurePath.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 1);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0}));

    closureParameters = GraphClosureParameters({.input_set = {0, 49}});
    graphClosurePath.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 50);
    closureParameters = GraphClosureParameters({.input_set = {1, 3}});
    graphClosurePath.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 3);
    EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({1, 2, 3}));

    GraphStruct circle = SimplePatterns::Circle(4);
    GraphClosure graphClosureCircle = GraphClosure(circle);
    closureParameters = GraphClosureParameters({.input_set = {0, 2}});
    graphClosureCircle.closure(closureParameters);
    EXPECT_EQ(closureParameters.closed_set.size(), 4);


    for (int i = 3; i < 100; ++i) {
        GraphStruct variableCircle = SimplePatterns::Circle(i);
        GraphClosure graphClosureVariableCircle = GraphClosure(variableCircle);
        closureParameters = GraphClosureParameters({.input_set = {0, (NodeId) i/2}});
        graphClosureVariableCircle.closure(closureParameters);
        // if i is odd
        if (i % 2 == 1) {
            EXPECT_EQ(closureParameters.closed_set.size(), (size_t) i/2 + 1);
        }
        else{
            EXPECT_EQ(closureParameters.closed_set.size(), i);
        }

        GraphStruct starGraph = SimplePatterns::StarGraph(i);
        GraphClosure graphClosureStarGraph = GraphClosure(starGraph);
        closureParameters = GraphClosureParameters({.input_set = {0, (NodeId) i - 1}});
        graphClosureStarGraph.closure(closureParameters);
        EXPECT_EQ(closureParameters.closed_set.size(), 2);
        EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, (NodeId) i-1}));
        closureParameters = GraphClosureParameters({.input_set = {1, (NodeId) i - 1}});
        graphClosureStarGraph.closure(closureParameters);
        EXPECT_EQ(closureParameters.closed_set.size(), 3);
        EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, 1, (NodeId) i-1}));
    }
}

TEST(GraphClosureTestSuite, ExampleClosureTree){

}

TEST(GraphClosureTestSuite, ExampleClosureOuterplanarGraph){

}


TEST(GraphClosureTestSuite, ExamplePreClosureTest){
    for (int i = 0; i < 20; ++i) {
        GraphStruct testGraph = SimplePatterns::MaxPreClosure(i);
        GraphClosure graphClosureTestGraph = GraphClosure(testGraph);
        for (int j = 0; j < testGraph.nodes(); ++j) {
            for (int k = 0; k < testGraph.nodes(); ++k) {
                GraphClosureParameters graphClosureParameters = GraphClosureParameters({.input_set = {(NodeId) j, (NodeId) k}});
                graphClosureTestGraph.closure(graphClosureParameters);
                if (j==k){
                    EXPECT_EQ(graphClosureParameters.closed_set.size(), 1);
                    EXPECT_EQ(graphClosureParameters.closed_set, std::set<NodeId>({(NodeId) j}));
                }
                else if (testGraph.edge(j, k)){
                    EXPECT_EQ(graphClosureParameters.closed_set.size(), 2);
                    EXPECT_EQ(graphClosureParameters.closed_set, std::set<NodeId>({(NodeId) j, (NodeId) k}));
                }
                else{
                    EXPECT_EQ(graphClosureParameters.closed_set.size(), testGraph.nodes());
                }
            }
        }
    }
}

TEST(GraphClosureTestSuite, ExampleClosureThetaTest){
    for (int i = 1; i <= 100; ++i) {
        GraphStruct testGraph = SimplePatterns::Path(i);
        GraphClosure graphClosureTestGraph = GraphClosure(testGraph);
        GraphClosureParameters closureParameters = GraphClosureParameters({.input_set = {0, (NodeId) i}, .threshold = 50});
        graphClosureTestGraph.closure(closureParameters);
        if (i <= 50){
            EXPECT_EQ(closureParameters.closed_set.size(), i + 1);
        }
        else{
            EXPECT_EQ(closureParameters.closed_set.size(), 2);
            EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({0, (NodeId) i}));
        }
    }

    GraphStruct testGraph = SimplePatterns::Path(9);
    GraphClosure graphClosureTestGraph = GraphClosure(testGraph);
    GraphClosureParameters closureParameters = GraphClosureParameters({.input_set = {1, 3, 6, 8}});
    for (int threshold = 0; threshold < 10; ++threshold) {
        closureParameters.threshold = threshold;
        graphClosureTestGraph.closure(closureParameters);
        if (threshold < 2){
            EXPECT_EQ(closureParameters.closed_set.size(), 4);
            EXPECT_EQ(closureParameters.closed_set, closureParameters.input_set);
        }
        else if (threshold == 2){
            EXPECT_EQ(closureParameters.closed_set.size(), 6);
            EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({1, 2, 3, 6, 7, 8}));
        }
        else {
            EXPECT_EQ(closureParameters.closed_set.size(), 8);
            EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({1, 2, 3, 4, 5, 6, 7, 8}));
        }
    }
    for (int threshold = 10; threshold >= 0; --threshold) {
        closureParameters.threshold = threshold;
        graphClosureTestGraph.closure(closureParameters);
        if (threshold < 2){
            EXPECT_EQ(closureParameters.closed_set.size(), 4);
            EXPECT_EQ(closureParameters.closed_set, closureParameters.input_set);
        }
        else if (threshold == 2){
            EXPECT_EQ(closureParameters.closed_set.size(), 6);
            EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({1, 2, 3, 6, 7, 8}));
        }
        else {
            EXPECT_EQ(closureParameters.closed_set.size(), 8);
            EXPECT_EQ(closureParameters.closed_set, std::set<NodeId>({1, 2, 3, 4, 5, 6, 7, 8}));
        }
    }

}

#endif //GOOGLE_TESTS_GRAPHCLOSURETESTS_H
