//
// Created by Florian on 16.03.2021.
//

#ifndef HOPS_EVALUATION_H
#define HOPS_EVALUATION_H

#include <iostream>
#include <utility>
#include <vector>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <map>
#include "HopsParameters.h"
#include "../GraphStaticFunctions.h"
#include "../io/StaticFunctions.h"
#include "../io/FileEvaluation.h"


class Evaluation {
public:
    Evaluation() = default;
    Evaluation(GraphStruct& graph, GraphStruct& pattern, RunParameters& runParameters) : graph(&graph), pattern(&pattern), parameters(runParameters), threads(runParameters.thread_num), seed(runParameters.seed){
        snapShotTimes = runParameters.runtime;
    };
    double HopsRuntime() const{return (double)(hopsRuntime.count())/1000000.0;};
    long double hopsEstimation = 0;
    long double hopsStd = 0;
    int threads = 0;
    int seed = 0;
    std::map<size_t, size_t> estimationMap;
    std::vector<double> snapShotTimes;
    std::vector<long double> snapshots;
    std::vector<long double> snapShotErrors;
    std::vector<UInt64> snapshotIterations;
    UInt64 hopsIterations = 0;
    UInt64 hopsZeroIterations = 0;

    //Time vars
    std::chrono::microseconds loadingTime{};
    std::chrono::microseconds preprocessingTime{};
    std::chrono::microseconds hopsRuntime{};
    std::chrono::microseconds evaluationRuntime{};

    //Run parameters
    RunParameters parameters;
private:
    GraphStruct* graph = nullptr;
    GraphStruct* pattern = nullptr;

    std::ofstream fs;




public:
    void print(HOPS_TYPE hopsType = HOPS_TYPE::ESTIMATION){
        switch (hopsType) {
            case HOPS_TYPE::ESTIMATION:
                std::cout << "////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
                std::cout << "\t" << "Graph: " << graph->GetName() << std::endl;
                std::cout << "\t" << "Pattern: " << pattern->GetName() << std::endl;
                if (parameters.labelType != LABEL_TYPE::UNLABELED) {
                    std::cout << "\t" << "Labeled embedding" << std::endl;
                }
                else{
                    std::cout << "\t" << "Unlabeled embedding" << std::endl;
                }
                std::cout << "\t" << "Estimation: " << hopsEstimation << std::endl;
                std::cout << "\t" << "HopsVariance: " << hopsStd << std::endl;
                std::cout << "\t" << "Iterations per Node: " << parameters.iteration_per_node.back() << std::endl;
                std::cout << "\t" << "Overall Iterations: " << hopsIterations << std::endl;
                std::cout << "\t" << "NonZeroIterations: " << hopsIterations - hopsZeroIterations << std::endl;
                std::cout << "\t" << "ZeroIterations: " << hopsZeroIterations << std::endl;
                std::cout << "---Parallelization---" << std::endl;
                std::cout << "\t" << "Threads: " << parameters.thread_num << std::endl;
                std::cout << "---Runtime Performance---" << std::endl;
                std::cout << "\t" << "HopsRuntime: " << (double) hopsRuntime.count() / 1000000 << " seconds" << std::endl;
                std::cout << "\t" << "Iterations/second: " << (double) hopsIterations / ((double) hopsRuntime.count() / 1000000) << std::endl;
                std::cout << "\t" << "Iterations/(Threads * second): " << (double) hopsIterations / ((double) hopsRuntime.count() / 1000000 * parameters.thread_num) << std::endl;
                std::cout << "\t" << "NonZeroIterations/second: " << (double) (hopsIterations - hopsZeroIterations) / ((double) hopsRuntime.count() / 1000000) << std::endl;
                std::cout << "\t" << "PreprocessingTime: " << (double) preprocessingTime.count() / 1000000 << " seconds" << std::endl;
                std::cout << "\t" << "EvaluationRuntime: " << (double) evaluationRuntime.count() / 1000000 << " seconds" << std::endl;

                std::cout << "---Additional Evaluation---" << std::endl;
                std::cout << "\t" << "Snapshot Times: " << std::endl;
                StaticFunctionsLib::print(snapShotTimes);
                std::cout << "\t" << "Snapshot Iterations: " << std::endl;
                StaticFunctionsLib::print(snapshotIterations);
                std::cout << "\t" << "Snapshot Values: " << std::endl;
                StaticFunctionsLib::print(snapshots);
                std::cout << "\t" << "Snapshots Errors: " << std::endl;
                for (size_t i = 1; i < snapshots.size(); ++i) {
                    snapShotErrors.emplace_back(std::abs(snapshots[i-1] - snapshots[i])/snapshots[i]);
                }
                StaticFunctionsLib::print(snapShotErrors);
                std::cout << "\t" << "Distinct Estimation Values: " << std::endl;
                StaticFunctionsLib::print(estimationMap);
                break;
            case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
                std::cout << "////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
                std::cout << "\t" << "Graph: " << graph->GetName() << std::endl;
                std::cout << "\t" << "Pattern: " << pattern->GetName() << std::endl;
                if (parameters.labelType != LABEL_TYPE::UNLABELED) {
                    std::cout << "\t" << "Labeled embedding" << std::endl;
                }
                else{
                    std::cout << "\t" << "Unlabeled embedding" << std::endl;
                }
                std::cout << "\t" << "ApproximatedGED: " << hopsEstimation << std::endl;
                std::cout << "\t" << "Iterations: " << hopsIterations << std::endl;
                std::cout << "---Parallelization---" << std::endl;
                std::cout << "\t" << "Threads:" << parameters.thread_num << std::endl;
                std::cout << "---Runtime Performance---" << std::endl;
                std::cout << "\t" << "HopsRuntime: " << (double) hopsRuntime.count() / 1000000 << " seconds" << std::endl;
                std::cout << "\t" << "Iterations/second: " << (double) hopsIterations / ((double) hopsRuntime.count() / 1000000) << std::endl;
                std::cout << "\t" << "Iterations/(Threads * second): " << (double) hopsIterations / ((double) hopsRuntime.count() / 1000000 * parameters.thread_num) << std::endl;
                std::cout << "\t" << "PreprocessingTime: " << (double) preprocessingTime.count() / 1000000 << " seconds" << std::endl;
                break;
        }

    };

    void save(const std::string& Path){
        FileEvaluation fileEvaluation = FileEvaluation(Path, this->graph->GetName(), ".hops");
        double hopsTime = (double) hopsRuntime.count() / 1000000;
        // get timestamp
        std::time_t t = std::time(nullptr);
        fileEvaluation.headerValueInsert({"TimeStamp", "Graph", "Size", "Edges", "PatternName", "PatternSize",
                                          "PatternEdges", "IsPatternTree", "Threads", "HopsTime",
                                          "HopsCount", "HopsIterationsPerNode", "HopsIterations", "Zero", "NonZero",
                                          "HopsIterations/second","Zero/second", "NonZero/second", "PreprocessingTime", "EvaluationTime", "SnapshotTimes", "SnapshotIterations",
                                          "Snapshots", "SnapshotErrors"},
                                         {std::to_string(t), graph->GetName(), std::to_string(graph->nodes()), std::to_string(graph->edges()), pattern->GetName(), std::to_string(pattern->nodes()),
                                          std::to_string(pattern->edges()), std::to_string(pattern->CheckTree()), std::to_string(threads), std::to_string(hopsTime),
                                          std::to_string(hopsEstimation), std::to_string(parameters.iteration_per_node.back()), std::to_string(hopsIterations), std::to_string(hopsZeroIterations), std::to_string(hopsIterations - hopsZeroIterations),
                                          std::to_string((double) hopsIterations/hopsTime), std::to_string((double ) hopsZeroIterations/hopsTime), std::to_string((double)(hopsIterations - hopsZeroIterations)/hopsTime),
                                          std::to_string((double) preprocessingTime.count() / 1000000), std::to_string(hopsTime ), StaticFunctionsLib::vectorToString(snapShotTimes), StaticFunctionsLib::vectorToString(snapshotIterations),
                                          StaticFunctionsLib::vectorToString(snapshots), StaticFunctionsLib::vectorToString(snapShotErrors)});
        fileEvaluation.save();
   };


};
#endif