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
        GraphClosure graphClosureTestGraph2 = GraphClosure(testGraph);
        for (int j = 0; j < testGraph.nodes(); ++j) {
            for (int k = 0; k < testGraph.nodes(); ++k) {
                GraphClosureParameters graphClosureParameters = GraphClosureParameters({.input_set = {(NodeId) j, (NodeId) k}});
                GraphClosureParameters graphClosureParameters2 = GraphClosureParameters({.input_set = {(NodeId) j, (NodeId) k}});
                graphClosureParameters2.closureType = EGraphClosureType::EXACT_GEODESIC_ITERATIVE;
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
        GraphClosureParameters closureParameters = GraphClosureParameters({.input_set = {0, (NodeId) i}, .theta = 50});
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
        closureParameters.theta = threshold;
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
        closureParameters.theta = threshold;
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

TEST(GraphClosureTestSuite, ExampleIterativeClosure)
{
    // generate a random graph with 10000 nodes and 1000000 edges
    GraphStruct randomGraph = SimplePatterns::ErdosRenyi(10000, 1000000);
    GraphClosure graphClosureRandomGraph = GraphClosure(randomGraph);

    std::set<NodeId> generator_set = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    // try different types of closures and report their running time
    GraphClosureParameters closureParameters = GraphClosureParameters({.input_set = generator_set});
    closureParameters.closureType = EGraphClosureType::EXACT_GEODESIC;

    // original closure
    auto start = std::chrono::high_resolution_clock::now();
    graphClosureRandomGraph.closure(closureParameters);
    auto end = std::chrono::high_resolution_clock::now();
    std::set<NodeId> original_closure = closureParameters.closed_set;
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Original Closure: " << elapsed_seconds.count() << "s\n";

    // iterative closure with runtime
    start = std::chrono::high_resolution_clock::now();
    closureParameters.closureType = EGraphClosureType::EXACT_GEODESIC_ITERATIVE;
    graphClosureRandomGraph.closure(closureParameters);
    end = std::chrono::high_resolution_clock::now();
    std::set<NodeId> iterative_closure = closureParameters.closed_set;
    elapsed_seconds = end-start;
    std::cout << "Iterative Closure: " << elapsed_seconds.count() << "s\n";

    // print sizes of closures
    std::cout << "Original Closure Size: " << original_closure.size() << "\n";
    std::cout << "Iterative Closure Size: " << iterative_closure.size() << "\n";

    EXPECT_EQ(original_closure, iterative_closure);
}

TEST(GraphClosureTestSuite, ExampleApproximateClosuresTrees)
{
    GraphExtended graph = GraphExtended(SimplePatterns::ErdosRenyi(1000, 20000, 0, true));
    GraphClosure graphClosure = GraphClosure(graph);
    std::vector<GraphStruct> subtrees;
    BFSSpanningTrees(graph, subtrees, 100, 0);
    GraphClosureParameters graphClosureParameters;
    graphClosureParameters.input_set = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    graphClosureParameters.closureType = EGraphClosureType::EXACT_GEODESIC;
    graphClosure.closure(graphClosureParameters);
    std::set<NodeId> exact_closure = graphClosureParameters.closed_set;

    graphClosureParameters.closureType = EGraphClosureType::APPROXIMATION_SUBSTRUCTURE;
    graphClosureParameters.clear();
    for (auto& x : subtrees){
        graphClosureParameters.sub_structures.emplace_back(&x);
    }
    graphClosure.closure(graphClosureParameters);

    std::set<NodeId> approx_closure = graphClosureParameters.closed_set;
    // print Jaccard Similarity between original closure and approximate closure
    NodeId intersection_size = 0;
    for (NodeId node : exact_closure)
    {
        if (approx_closure.find(node) != approx_closure.end())
        {
            intersection_size++;
        }
    }
    double jaccard_similarity = (double) intersection_size / (exact_closure.size() + approx_closure.size() - intersection_size);
    std::cout << "Jaccard Similarity: " << jaccard_similarity << "\n";
}

TEST(GraphClosureTestSuite, ExampleApproximateClosuresOuterplanarSubgraphs)
{
    GraphExtended graph = GraphExtended(SimplePatterns::ErdosRenyi(1000, 20000, 0, true));
    GraphClosure graphClosure = GraphClosure(graph);
    std::vector<GraphStruct> outerplanar_subgraphs;
    std::vector<OuterplanarGraphData> outerplanar_graph_data;
    OuterplanarSubgraphDFS outerplanarSubgraphDFS = OuterplanarSubgraphDFS(graph);
    outerplanarSubgraphDFS.subgraphs(100, outerplanar_subgraphs, 0);
    outerplanarSubgraphDFS.subgraphs_extended(100, outerplanar_graph_data, 0);
    GraphClosureParameters graphClosureParameters;
    graphClosureParameters.input_set = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    graphClosureParameters.closureType = EGraphClosureType::EXACT_GEODESIC;
    graphClosure.closure(graphClosureParameters);
    std::set<NodeId> exact_closure = graphClosureParameters.closed_set;

    graphClosureParameters.closureType = EGraphClosureType::APPROXIMATION_SUBSTRUCTURE;
    // fill substructures with pointers to outerplanar_subgraphs
    graphClosureParameters.sub_structures.clear();
    for (auto& x : outerplanar_subgraphs){
        graphClosureParameters.sub_structures.emplace_back(&x);
    }
    graphClosure.closure(graphClosureParameters);
    std::set<NodeId> approx_closure1 = graphClosureParameters.closed_set;

    graphClosureParameters.sub_structures.clear();
    for (auto& x : outerplanar_graph_data){
        graphClosureParameters.sub_structures.emplace_back(&x);
    }
    graphClosure.closure(graphClosureParameters);
    std::set<NodeId> approx_closure2 = graphClosureParameters.closed_set;

    EXPECT_NE(approx_closure1.size(), 0);
    EXPECT_EQ(approx_closure1, approx_closure2);

    // print Jaccard Similarity between original closure and approximate closure
    NodeId intersection_size = 0;
    for (NodeId node : exact_closure)
    {
        if (approx_closure1.find(node) != approx_closure1.end())
        {
            intersection_size++;
        }
    }
    double jaccard_similarity = (double) intersection_size / (exact_closure.size() + approx_closure1.size() - intersection_size);
    std::cout << "Jaccard Similarity: " << jaccard_similarity << "\n";
}

TEST(GraphClosureTestSuite, ExampleApproximateClosuresOuterplanarSubgraphsComparison)
{
    // load CA-GrQc graph
    std::string path = "../../../../GraphData/RealWorld/Collaboration/CA-GrQc_component.bin";
    // print files in path
    // check if file exists
    if (!std::filesystem::exists(path))
    {
        std::cout << "File does not exist\n";
        return;
    }

    GraphStruct graph = GraphStruct(path);
    OuterplanarSubgraphDFS outerplanarSubgraphDFS = OuterplanarSubgraphDFS(graph);
    GraphStruct subgraph;
    outerplanarSubgraphDFS.generate(subgraph, 0, false);
    OuterplanarGraphData outerplanarGraphData;
    outerplanarSubgraphDFS.subgraph_extended(subgraph, outerplanarGraphData, 0, false);

    GraphClosureParameters graphClosureParameters;
    // get 5 random nodes
    StaticFunctionsLib::get_k_from_n(graphClosureParameters.input_set, 5, graph.nodes(), 0);

    GraphClosure graphClosure = GraphClosure(subgraph);
    graphClosure.closure(graphClosureParameters);
    std::set<NodeId> closure = graphClosureParameters.closed_set;

    GraphClosure graphClosure2 = GraphClosure(outerplanarGraphData);
    graphClosure2.closure(graphClosureParameters);
    std::set<NodeId> closure2 = graphClosureParameters.closed_set;

    EXPECT_EQ(closure, closure2);

}


#endif //GOOGLE_TESTS_GRAPHCLOSURETESTS_H
