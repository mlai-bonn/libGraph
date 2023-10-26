//
// Created by florian on 26.10.23.
//

#ifndef GOOGLE_TESTS_GRAPHCORETESTS_H
#define GOOGLE_TESTS_GRAPHCORETESTS_H
#include <gtest/gtest.h>
#include "libGraph.h"

TEST(GraphCoreTestSuite, ExampleCoreRandomAlgorithmErdosRenyi){
    GraphStruct graph = SimplePatterns::ErdosRenyi(1000, 20000, 0, true);
    auto coreRandomAlgorithm = CoreRandomAlgorithm(graph);
    CoreAlgorithmParameters parameters = {.generator_size = 5, .output_path = "../tests/out/"};
    coreRandomAlgorithm.Run(parameters);
    EXPECT_EQ(parameters.core_nodes.size(), graph.nodes());

    parameters = {.generator_size = 5, .core_iterations = 10, .output_path = "../tests/out/"};
    coreRandomAlgorithm.Run(parameters);
    EXPECT_EQ(parameters.core_nodes.size(), graph.nodes());

    parameters = {.generator_size = 1};
    coreRandomAlgorithm.Run(parameters);
    EXPECT_EQ(parameters.core_nodes.size(), 0);
}

TEST(GraphCoreTestSuite, ExampleCoreGrowAlgorithmErdosRenyi){
    GraphStruct graph = SimplePatterns::ErdosRenyi(1000, 20000, 0, true);
    auto coreGrowAlgorithm = CoreGrowAlgorithm(graph);
    CoreGrowAlgorithmParameters parameters = {.num_runs = 10, .grow_steps = 5, .core_percentage = 0.9, .output_path = "../tests/out/"};
    coreGrowAlgorithm.Run(parameters);
    EXPECT_EQ(parameters.core_nodes.size(), graph.nodes());
}



#endif //GOOGLE_TESTS_GRAPHCORETESTS_H
