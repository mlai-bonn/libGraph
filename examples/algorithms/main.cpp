// This file shows how to use the libGraph graph library


#include "SimplePatterns.h"
#include "test.h"
#include "Algorithms/GED/GEDApproximation.h"
#include "GraphDataStructures/GraphBase.h"

// Example for calculating an approximation of the Graph Edit Distance between two graphs
GEDApproximationParameters graph_edit_distance_approximation() {
    GraphStruct start_graph = SimplePatterns::Circle(3);
    GraphStruct target_graph = SimplePatterns::Circle(4);



    GEDApproximationParameters parameters = {
        .type = GEDApproximationType::ASTAR,
        .seed = 0,

    };

    GEDApproximation ged_approximation = GEDApproximation(start_graph, target_graph);
    ged_approximation.Run(parameters);
    return parameters;
}

// main function to run the examples

int main() {
    test();
    GEDApproximationParameters parameters = graph_edit_distance_approximation();
    return 0;
}

