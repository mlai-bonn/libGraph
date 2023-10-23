//
// Created by florian on 17.08.23.
//

#ifndef TESTGRAPHLIB_COREALGORITHM_H
#define TESTGRAPHLIB_COREALGORITHM_H

#include "../../typedefs.h"
#include "../../GraphDataStructures/GraphBase.h"

struct CoreAlgorithmInputParameters {
    int generator_size = 5;
    int core_iterations = -1;
    int seed = 0;
    bool print = false;
    bool save = false;
};

struct CoreAlgorithmOutputParameters {
    std::vector<NodeId> core_nodes = std::vector<NodeId>();
    std::vector<int> intersection_loss = std::vector<int>();
    double runtime = 0.0;
};

class CoreAlgorithm {
public:
    explicit CoreAlgorithm(GraphStruct &graph) : _graph(graph) {};
    void Run(CoreAlgorithmOutputParameters& outputParameters, const CoreAlgorithmInputParameters& inputParameters){
        GraphClosureSP gc = GraphClosureSP(_graph);
        ClosureParameters closureParameters;
        std::vector<std::set<NodeId>> closures;
        std::cout << std::endl;
        std::set<NodeId> overlap;
        outputParameters.runtime = 0;
        int c_iterations = inputParameters.core_iterations;
        if (inputParameters.core_iterations == -1) {
            c_iterations = _graph.nodes();
        }

        // start timer
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < c_iterations; ++i) {
            std::cout << "\tIteration " << std::to_string(i) << " of graph " << _graph.GetName() << " started"
                      << std::endl;

            // generate input set
            std::mt19937_64 generator(i + c_iterations * inputParameters.seed);
            closureParameters.input_set.clear();
            std::vector<NodeId> nodes(_graph.nodes());
            std::iota(nodes.begin(), nodes.end(), 0);
            for (int i = 0; i < inputParameters.generator_size; ++i) {
                int rand_idx = std::uniform_int_distribution<int>(i, ((int) nodes.size()) - 1)(generator);
                closureParameters.input_set.insert(nodes[rand_idx]);
                std::swap(nodes[rand_idx], nodes[i]);
            }

            // compute the closure
            gc.closure(closureParameters);
            closures.emplace_back(closureParameters.closed_set);
            if (i == 0) {
                overlap = closureParameters.closed_set;
                std::cout << "\tIteration " << std::to_string(i) << " of graph " << _graph.GetName() << " finished after "
                          << std::to_string(((double) std::chrono::duration_cast<std::chrono::microseconds>(
                                  std::chrono::high_resolution_clock::now() - start).count() / 1000000.0)) << "s"
                          << std::endl << std::endl;
            } else {
                int overlap_size = (int) overlap.size();
                std::vector<NodeId> v_intersection;
                std::set_intersection(overlap.begin(), overlap.end(),
                                      closures.back().begin(), closures.back().end(),
                                      std::back_inserter(v_intersection));
                overlap.clear();
                overlap.insert(v_intersection.begin(), v_intersection.end());
                outputParameters.intersection_loss.emplace_back(overlap_size - overlap.size());
                std::cout << "\tIteration " << std::to_string(i) << " of graph " << _graph.GetName() << " finished after "
                          << std::to_string(((double) std::chrono::duration_cast<std::chrono::microseconds>(
                                  std::chrono::high_resolution_clock::now() - start).count() / 1000000.0)) << "s"
                          << std::endl << std::endl;
                if (inputParameters.core_iterations == -1 && outputParameters.intersection_loss.back() == 0) {
                    break;
                }
            }

        }
        outputParameters.core_nodes.clear();
        for (auto elem: overlap) {
            outputParameters.core_nodes.emplace_back(elem);
        }
        //Measure runtime before getting statistics
        outputParameters.runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count() /
                                    1000000.0);
    }

    GraphStruct& _graph;



};


#endif //TESTGRAPHLIB_COREALGORITHM_H
