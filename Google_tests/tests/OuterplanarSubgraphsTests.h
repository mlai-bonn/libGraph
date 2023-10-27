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
    GraphStruct circle_diagonals = SimplePatterns::Circle(1000);
    for(int i = 1; i < 999; ++i){
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
    OuterPlanarSubgraphMitchell outerPlanarSubgraphMitchell = OuterPlanarSubgraphMitchell(fullyBipartite);
    GraphStruct o1 = outerplanarSubgraphDFS.subgraph(0, false);
    GraphStruct o2 = outerPlanarSubgraphMitchell.subgraph(0, false);
    EXPECT_EQ(o1.edges(), 5);
    EXPECT_EQ(o2.edges(), 5);
}

#endif //GOOGLE_TESTS_OUTERPLANARSUBGRAPHSTESTS_H
