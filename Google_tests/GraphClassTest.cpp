//
// Created by florian on 24.10.23.
//

#include <gtest/gtest.h>
#include "GraphDataStructures/GraphBase.h"

TEST(GraphConstructorTestSuite, ExampleGraphConstructor){
GraphStruct graph = GraphStruct();

EXPECT_EQ(graph.nodes(), 0);
}