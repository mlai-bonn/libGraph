//
// Created by florian on 17.10.23.
//
// This file contains the algorithm from the paper "Constant approximation algorithms for embedding
//graph metrics into trees and outerplanar graphs" by Chepoi, Dragan, Newman, Rabinovich and Vaxès.

#ifndef LIBGRAPH_TREEEMBEDDING_H
#define LIBGRAPH_TREEEMBEDDING_H

#include "Layering.h"
#include "../Graph/GraphAlgorithms.h"
#include "Algorithms/Graph/SpanningTrees.h"

enum class ApproximationType{
    RANDOM_SPANNING_TREES,
    OUTERPLANAR_SPANNING_GRAPHS,
    LAYERING_TREE,
};

struct ParametersEvaluation{
    ApproximationType approximationType = ApproximationType::LAYERING_TREE;
    int seed = 0;
    int root_node_id = -1;
    int tree_number = 1;
    bool store_all = false;
};



class EvaluateApproximations { ;

public:
    explicit EvaluateApproximations(GraphStruct &graph) : _graph(graph) {}

    void Evaluation(ParametersEvaluation& params);

private:
    GraphStruct& _graph;
};

void EvaluateApproximations::Evaluation(ParametersEvaluation &params) {
    // for each node run a BFS in the graph and in the layering trees and compare the results
    std::mt19937_64 generator(params.seed);
    std::vector<INDEX> graphDistances;
    std::vector<INDEX> approximationDistances;
    LayeringTree layeringTree;
    std::vector<GraphStruct> spanningTrees;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    // make a case distinction for the approximation type
    switch (params.approximationType) {
        case ApproximationType::RANDOM_SPANNING_TREES:
            // run the random spanning tree approximation
            for (int i = 0; i < params.tree_number; ++i){
                spanningTrees.emplace_back();
                NodeId root_node = std::uniform_int_distribution<NodeId>(0, _graph.nodes() - 1)(generator);
                BFSSpanningTree(_graph, spanningTrees.back(), root_node, visited, distances, false, params.seed + i);
            }
            break;
        case ApproximationType::OUTERPLANAR_SPANNING_GRAPHS:
            // run the outerplanar spanning graph approximation
            break;
        case ApproximationType::LAYERING_TREE:
            // run the layering tree approximation
            // create the layering trees
            layeringTree = LayeringTree(_graph, params.seed, params.root_node_id);
            break;
    }

    std::vector<std::vector<std::pair<int, int>>> distance_matrix;
    int false_distances = 0;
    int correct_distances = 0;

    double relative_error = 0.0;
    size_t absolute_error = 0;

    double mean_relative_error = 0.0;
    double mean_error = 0.0;

    int minimum_error = std::numeric_limits<int>::max();
    int maximum_error = std::numeric_limits<int>::min();


    if (params.store_all){
        distance_matrix = std::vector<std::vector<std::pair<int, int>>>(_graph.nodes(), std::vector<std::pair<int, int>>(_graph.nodes(), {0, 0}));
    }

    for (int i = 0; i < _graph.nodes(); ++i) {
        GraphStruct::BFSDistances(_graph, i, graphDistances);
        switch (params.approximationType) {
            case ApproximationType::RANDOM_SPANNING_TREES:
                std::fill(approximationDistances.begin(), approximationDistances.end(), std::numeric_limits<INDEX>::max());
                approximationDistances.resize(_graph.nodes(), std::numeric_limits<INDEX>::max());

                for (auto& spanningTree : spanningTrees){
                    std::vector<INDEX> spanningTreeDistances;
                    GraphStruct::BFSDistances(spanningTree, i, spanningTreeDistances);
                    for (int j = 0; j < _graph.nodes(); ++j) {
                        approximationDistances[j] = std::min(approximationDistances[j], spanningTreeDistances[j]);
                    }
                }
                // run the random spanning tree approximation
                break;
            case ApproximationType::OUTERPLANAR_SPANNING_GRAPHS:
                // run the outerplanar spanning graph approximation
                break;
            case ApproximationType::LAYERING_TREE:
                // run the layering tree approximation
                // create the layering trees
                GraphStruct::BFSDistances(layeringTree.tree, i, approximationDistances);
                break;
        }

        for (int j = 0; j < _graph.nodes(); ++j) {
            // count only non-zero distances
            if (graphDistances[j] != 0) {
                if (graphDistances[j] != approximationDistances[j]) {
                    ++false_distances;
                } else {
                    ++correct_distances;
                }
                if (params.store_all) {
                    distance_matrix[i][j] = {graphDistances[j], approximationDistances[j]};
                }
                if (graphDistances[j] > approximationDistances[j]) {
                    absolute_error += graphDistances[j] - approximationDistances[j];
                } else {
                    absolute_error += approximationDistances[j] - graphDistances[j];

                }
                minimum_error = std::min(minimum_error, (int) graphDistances[j] - (int) approximationDistances[j]);
                maximum_error = std::max(maximum_error, (int) graphDistances[j] - (int) approximationDistances[j]);
                relative_error += (double) approximationDistances[j] / (double) graphDistances[j];
            }
        }
    }


    mean_error = (double) absolute_error / (false_distances + correct_distances);
    mean_relative_error = relative_error / (false_distances + correct_distances);

    // print the results
    // print the root node
    switch (params.approximationType) {

        case ApproximationType::RANDOM_SPANNING_TREES:
            std::cout << "Random Spanning Tree Evaluation" << std::endl;
            std::cout << "Number of spanning trees: " << params.tree_number << std::endl;
            break;
        case ApproximationType::OUTERPLANAR_SPANNING_GRAPHS:
            break;
        case ApproximationType::LAYERING_TREE:
            std::cout << "Root node: " << layeringTree.layeringPartition.root_node_id() << std::endl;
            std::cout << "Layering Tree Evaluation" << std::endl;
            break;
    }

    std::cout << "Number of false distances: " << false_distances << std::endl;
    std::cout << "Number of correct distances: " << correct_distances << std::endl;
    std::cout << "Absolute error: " << absolute_error << std::endl;
    std::cout << "Mean absolute error: " << mean_error << std::endl;

    std::cout << "Mean relative error: " << mean_relative_error << std::endl;

    std::cout << "Minimum error: " << minimum_error << std::endl;
    std::cout << "Maximum error: " << maximum_error << std::endl;

    // store all values to file using FileEvaluation
    FileEvaluation evaluation = FileEvaluation();
    std::string filename;
    switch (params.approximationType) {
        case ApproximationType::RANDOM_SPANNING_TREES:
            // run the random spanning tree approximation
            evaluation = FileEvaluation("../out/", "random_spanning_trees_evaluation");
            evaluation.headerValueInsert({"Number of Spanning Trees", "False Distances", "Correct Distances", "Absolute Error", "Mean Absolute Error", "Mean Relative Error", "Minimum Error", "Maximum Error"},
                                         {std::to_string(params.tree_number), std::to_string(false_distances), std::to_string(correct_distances), std::to_string(absolute_error), std::to_string(mean_error), std::to_string(mean_relative_error), std::to_string(minimum_error), std::to_string(maximum_error)});
            filename = "../out/random_spanning_trees_distance_matrix.csv";
            break;
        case ApproximationType::OUTERPLANAR_SPANNING_GRAPHS:
            // run the outerplanar spanning graph approximation
            break;
        case ApproximationType::LAYERING_TREE:
            // run the layering tree approximation
            evaluation = FileEvaluation("../out/", "layering_tree_evaluation");
            evaluation.headerValueInsert({"Root Node", "False Distances", "Correct Distances", "Absolute Error", "Mean Absolute Error", "Mean Relative Error", "Minimum Error", "Maximum Error"},
                                         {std::to_string(layeringTree.layeringPartition.root_node_id()), std::to_string(false_distances), std::to_string(correct_distances), std::to_string(absolute_error), std::to_string(mean_error), std::to_string(mean_relative_error), std::to_string(minimum_error), std::to_string(maximum_error)});

            filename = "../out/layering_tree_distance_matrix.csv";
            // create the layering trees
            break;
    }

    evaluation.save();


    // if store all is true save error matrix to csv file in the out folder
    if (params.store_all){

        std::ofstream file;
        file.open(filename);
        for (int i = 0; i < _graph.nodes(); ++i){
            for (int j = 0; j < _graph.nodes(); ++j){
                file << "(" << distance_matrix[i][j].first << ";" << distance_matrix[i][j].second << ")" << ",";
            }
            file << std::endl;
        }
        file.close();
    }
}





#endif //LIBGRAPH_TREEEMBEDDING_H
