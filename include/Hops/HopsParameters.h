//
// Created by florian on 27.01.22.
//

#ifndef HOPS_HOPSPARAMETERS_H
#define HOPS_HOPSPARAMETERS_H

#include "string"
#include "vector"
#include "set"

enum class TestType{
    NONE,
    TRIANGLE_UNLABELED,
    TRIANGLE_LABELED_SPARSE,
    TRIANGLE_LABELED_DENSE
};

struct Parameters{
    /// parameters for the hops class, mainly the paths for loading and saving
    std::string in_path = "../../GraphData/";
    std::string out_path = "../results/results.csv";
    std::vector<std::string> graph_names = {"amazon"};
    std::string pattern_path = "../data/patterns/" + graph_names[0] + "/cyclic_patterns/";
    TestType testType = TestType::NONE;
    std::string test_mode;
};

struct RunParameters{
    /// parameters for one run (one pattern count estimation) of the hops algorithm
    LABEL_TYPE labelType = LABEL_TYPE::UNLABELED;
    double runtime = 30;
    int thread_num = -1;
    std::vector<int> iteration_per_node = {0};
    int seed = 0;
    bool print = true;
    bool save = false;
    bool single_number = true;
    int snapshot_time = -1;

    HOPS_TYPE hops_type = HOPS_TYPE::ESTIMATION;
    RUN_TYPE run_type = RUN_TYPE::VECTOR;


    DGraphStruct spanningTree = DGraphStruct();
    NodeId spanningTreeRootNode = std::numeric_limits<NodeId>::max();
};


#endif //HOPS_HOPSPARAMETERS_H
