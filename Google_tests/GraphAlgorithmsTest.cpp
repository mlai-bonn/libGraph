// Testing the graph algorithms

#include <gtest/gtest.h>
#include "../include/libGraph.h"

TEST(BiconnectedComponentTestSuite, ExampleBiconnectedComponent){
    GraphStruct outerplanar_graph = Chepoi_Outerplanar_Graph();
    std::vector<std::vector<NodeId>> components;
    GetBiconnectedComponents(outerplanar_graph, components);

    EXPECT_EQ(components.size(), 1);
    
    GraphStruct bi_conn_wiki_graph = Bi_Conn_Wiki_Graph();

    GetBiconnectedComponents(bi_conn_wiki_graph, components);
    EXPECT_EQ(components.size(), 7);

}

TEST(OuterplanarGraphDataTestSuite, ExampleOuterplanarGraphData){
    GraphStruct bi_conn_wiki_graph = Bi_Conn_Wiki_Graph();
    OuterplanarGraphData outerplanarGraphData = OuterplanarGraphData(bi_conn_wiki_graph, 14);
    outerplanarGraphData.set();
}

