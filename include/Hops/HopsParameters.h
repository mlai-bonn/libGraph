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
    /// parameters for a run of the hops algorithm
    LABEL_TYPE labelType = LABEL_TYPE::UNLABELED;
    double runtime = 30;
    int thread_num = -1;
    int iteration_per_thread = 0;
    int seed = 0;
    bool print = true;
    bool single_number = true;
    int snapshot_time = -1;

    HOPS_TYPE hops_type = HOPS_TYPE::ESTIMATION;
    RUN_TYPE run_type = RUN_TYPE::VECTOR;

    std::string run_type_string;

    std::set<int> pattern_sizes = {4};
    int pattern_num = 1;

    DGraphStruct spanningTree = DGraphStruct();
    INDEX spanningTreeRootNode = std::numeric_limits<INDEX>::max();


};


#endif //HOPS_HOPSPARAMETERS_H
