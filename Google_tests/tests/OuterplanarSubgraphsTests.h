//
// Created by florian on 26.10.23.
//

#ifndef GOOGLE_TESTS_OUTERPLANARSUBGRAPHSTESTS_H
#define GOOGLE_TESTS_OUTERPLANARSUBGRAPHSTESTS_H

TEST(OuterplanarSubgraphsTestSuite, ExampleOuterplanarDFSCircle){
    GraphStruct circle = SimplePatterns::Circle(4);

    OuterplanarSubgraphDFS outerplanarSubgraphDFS = OuterplanarSubgraphDFS(circle);
    OuterPlanarSubgraphMitchell outerPlanarSubgraphMitchell = OuterPlanarSubgraphMitchell(circle);
    GraphStruct& o1 = outerplanarSubgraphDFS.subgraph(0, false);
    //GraphStruct& o2 = outerPlanarSubgraphMitchell.subgraph(0, false);
    EXPECT_EQ(o1.nodes(), circle.nodes());
    EXPECT_EQ(o1.edges(), circle.edges());
    //EXPECT_EQ(o2.nodes(), circle.nodes());
    //EXPECT_EQ(o2.edges(), circle.edges());
    // the same edges
    for(auto edge = circle.first_edge(); edge != circle.last_edge(); ++edge){
        EXPECT_EQ(o1.edge((*edge).first, (*edge).second), true);
        //EXPECT_EQ(o2.edge((*edge).first, (*edge).second), true);
    }
}

TEST(OuterplanarSubgraphsTestSuite, ExampleOuterplanarDFSCircleWithDiagonals){
    GraphStruct circle_diagonals = SimplePatterns::Circle(10);
    for(int i = 1; i < 9; ++i){
        circle_diagonals.add_edge(0, i);
    }

    OuterplanarSubgraphDFS outerplanarSubgraphDFS = OuterplanarSubgraphDFS(circle_diagonals);
    GraphStruct& o1 = outerplanarSubgraphDFS.subgraph(0, false);
    //GraphStruct& o2 = outerPlanarSubgraphMitchell.subgraph(0, false);
    EXPECT_EQ(o1.nodes(), circle_diagonals.nodes());
    EXPECT_EQ(o1.edges(), circle_diagonals.edges());
    //EXPECT_EQ(o2.nodes(), circle.nodes());
    //EXPECT_EQ(o2.edges(), circle.edges());
    // the same edges
    for(auto edge = circle_diagonals.first_edge(); edge != circle_diagonals.last_edge(); ++edge){
        EXPECT_EQ(o1.edge((*edge).first, (*edge).second), true);
        //EXPECT_EQ(o2.edge((*edge).first, (*edge).second), true);
    }
}


TEST(OuterplanarSubgraphsTestSuite, ExampleOuterplanarDFSDoubleTriangle){
    GraphStruct doubleTriangle = SimplePatterns::DoubleTriangle();

    OuterplanarSubgraphDFS outerplanarSubgraphDFS = OuterplanarSubgraphDFS(doubleTriangle);
    OuterPlanarSubgraphMitchell outerPlanarSubgraphMitchell = OuterPlanarSubgraphMitchell(doubleTriangle);
    GraphStruct& o1 = outerplanarSubgraphDFS.subgraph(0, false);
    //GraphStruct& o2 = outerPlanarSubgraphMitchell.subgraph(0, false);
    EXPECT_EQ(o1.nodes(), doubleTriangle.nodes());
    EXPECT_EQ(o1.edges(), doubleTriangle.edges());
    //EXPECT_EQ(o2.nodes(), doubleTriangle.nodes());
    //EXPECT_EQ(o2.edges(), doubleTriangle.edges());
    // the same edges
    for(auto edge = doubleTriangle.first_edge(); edge != doubleTriangle.last_edge(); ++edge){
        EXPECT_EQ(o1.edge((*edge).first, (*edge).second), true);
        //EXPECT_EQ(o2.edge((*edge).first, (*edge).second), true);
    }
}

TEST(OuterplanarSubgraphsTestSuite, ExampleOuterplanarDFSFullyBipartite){
    GraphStruct fullyBipartite = SimplePatterns::FullyBipartite(2,3);
    OuterplanarSubgraphDFS outerplanarSubgraphDFS = OuterplanarSubgraphDFS(fullyBipartite);
    //OuterPlanarSubgraphMitchell outerPlanarSubgraphMitchell = OuterPlanarSubgraphMitchell(fullyBipartite);
    GraphStruct o1 = outerplanarSubgraphDFS.subgraph(0, false);
    //GraphStruct o2 = outerPlanarSubgraphMitchell.subgraph(0, false);
    EXPECT_EQ(o1.edges(), 5);
    //EXPECT_EQ(o2.edges(), 5);
}

TEST(OuterplanarSubgraphDataTestSuite, ExampleOuterplanarSubgraphData){
    int size = 3;
    int circle_number = 4;
    GraphStruct n_circles = GraphStruct(circle_number*(size - 1) + 1, {});
    for(int i = 0; i < circle_number; ++i){
        for(int j = 0; j < size; ++j){
            n_circles.add_edge(i*(size - 1) + j, i*(size - 1) + (j + 1) % size);
        }
    }
    OuterplanarGraphData outerplanarGraphData = OuterplanarGraphData(n_circles);
    EXPECT_EQ(outerplanarGraphData.nodes(), n_circles.nodes());
    EXPECT_EQ(outerplanarGraphData.edges(), n_circles.edges());
    EXPECT_EQ(outerplanarGraphData.Components.size(), circle_number);
}

TEST(OuterplanarSubgraphDataTestSuite, ExampleDetailedTest) {
    int seed = 0;
    int nodes = 10;
    int edges = 15;
    std::vector<int> neighbors = std::vector<int>();
    int nonOuterPlanarNumber = 0;
    std::vector<int> nonOuterPlanarSeeds;
    int nonMaximalNumber = 0;
    std::vector<int> nonMaximalSeeds;

    std::vector<double> maximalEdges;
    std::vector<double> nonMaximalEdgeNumPerCent;
    std::vector<double> outerPlanarEdgeNum;
    std::vector<double> maximalEdgesNum;
    std::vector<double> algorithmMissingEdges;
    std::vector<double> missingEdgesNum;
    size_t outerplanar_runtime = 0;
    size_t outerplanar_structure_runtime = 0;
    OuterplanarGraphStatistics statistics = OuterplanarGraphStatistics();
    auto graph = SimplePatterns::ErdosRenyi(nodes, edges, seed, true);
    // save graph as edge list to file
    graph.Save({.graphPath = "../tests/out/", .Name = "erdos_renyi", .Format = GraphFormat::EDGES});

    auto start = std::chrono::high_resolution_clock::now();
    auto subgraphGeneration = OuterplanarSubgraphDFS(graph);
    GraphStruct o1 = subgraphGeneration.subgraph(seed, false);

    outerplanar_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count();

    start = std::chrono::high_resolution_clock::now();
    OuterplanarGraphData outerplanarSubgraph = OuterplanarGraphData(o1);
    outerplanar_structure_runtime += std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start).count();
    outerplanarSubgraph.Save({.graphPath = "../tests/out/", .Name = "erdos_renyi_outerplanar_subgraph", .Format = GraphFormat::EDGES});

    //Check outerplanarity + quality
    GraphAlgorithms::CheckingOuterpanarity(graph, o1, nonOuterPlanarNumber,
                          nonMaximalNumber,
                          nonOuterPlanarSeeds,
                          nonMaximalSeeds, algorithmMissingEdges, maximalEdges,
                          seed);
    maximalEdgesNum.emplace_back(o1.edges() + algorithmMissingEdges.back());
    nonMaximalEdgeNumPerCent.emplace_back(
            (double) (o1.edges()) /
            ((double) o1.edges() + algorithmMissingEdges.back()) * 100.0);
    missingEdgesNum.emplace_back(algorithmMissingEdges.back());

    statistics += OuterplanarGraphStatistics(o1);
    outerPlanarEdgeNum.emplace_back(o1.edges());
    EXPECT_EQ(GraphStruct::IsConnected(graph), true);
    EXPECT_EQ(nonOuterPlanarNumber, 0);
}

#endif //GOOGLE_TESTS_OUTERPLANARSUBGRAPHSTESTS_H
