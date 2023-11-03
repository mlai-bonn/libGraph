//
// Created by florian on 17.08.23.
//

#ifndef TESTGRAPHLIB_COREALGORITHM_H
#define TESTGRAPHLIB_COREALGORITHM_H

/**
 * @brief Input parameters for the core algorithm introduced by Seiffarth et al. in their paper "A Fast Heuristic for Computing Geodesic Cores in Large Networks" (2022)
 * The algorithm is a probabilistic algorithm that finds the geodesic core of a graph, in each run it chooses some random nodes and then calculates the closure of the set of nodes. The algorithm stops either if the intersection of the closed sets computed so far gets stable (core iterations = -1) or if the number of iterations is reached (core iterations = n).
 * generator_size: size of the random input set
 * core_iterations: number of iterations of the algorithm, if -1 the algorithm stops if the intersection of the closed sets computed so far gets stable
 * _seed: _seed for the random number generator
 * _print: if true prints the growth steps and the closure size
 * save: if true saves the core nodes vector
 * output_path: path to the output file
 * core_nodes: vector of core nodes
 * intersection_loss: vector of the intersection loss in each iteration
 * runtime: runtime of the algorithm
 */
struct CoreAlgorithmParameters {

    // Input parameters
    int generator_size = 5;
    int core_iterations = -1;
    int seed = 0;
    bool print = false;
    bool save = true;


    // Output parameters
    FileEvaluation core_evaluation;
    FileEvaluation detailed_evaluation;
    std::vector<NodeId> core_nodes = std::vector<NodeId>();
    std::vector<int> intersection_loss = std::vector<int>();
    double runtime = 0.0;
};

class CoreRandomAlgorithm {
public:
    explicit CoreRandomAlgorithm(GraphStruct &graph) : _graph(graph) {};
    void Run(CoreAlgorithmParameters& parameters){
        GraphClosure gc = GraphClosure(_graph);
        GraphClosureParameters closureParameters;
        std::vector<std::set<NodeId>> closures;
        std::cout << std::endl;
        std::set<NodeId> overlap;
        parameters.runtime = 0;
        INDEX c_iterations = parameters.core_iterations;
        if (parameters.core_iterations == -1) {
            c_iterations = _graph.nodes();
        }

        // save the results
        parameters.core_evaluation = FileEvaluation();
        parameters.detailed_evaluation = FileEvaluation();

        // save the details of the algorithm
        // start timer
        auto start = std::chrono::high_resolution_clock::now();
        int iteration_count = 0;
        for (int i = 0; i < c_iterations; ++i) {
            std::cout << "\tIteration " << std::to_string(i) << " of _graph " << _graph.GetName() << " started"
                      << std::endl;

            // generate input set
            std::mt19937_64 generator(i + c_iterations * parameters.seed);
            closureParameters.input_set.clear();
            std::vector<NodeId> nodes(_graph.nodes());
            std::iota(nodes.begin(), nodes.end(), 0);
            for (int j = 0; j < parameters.generator_size; ++j) {
                int rand_idx = std::uniform_int_distribution<int>(j, ((int) nodes.size()) - 1)(generator);
                closureParameters.input_set.insert(nodes[rand_idx]);
                std::swap(nodes[rand_idx], nodes[j]);
            }

            // compute the closure
            gc.closure(closureParameters);
            closures.emplace_back(closureParameters.closed_set);
            if (i == 0) {
                overlap = closureParameters.closed_set;
                std::cout << "\tIteration " << std::to_string(i) << " of _graph " << _graph.GetName() << " finished after "
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
                parameters.intersection_loss.emplace_back(overlap_size - overlap.size());
                std::cout << "\tIteration " << std::to_string(i) << " of _graph " << _graph.GetName() << " finished after "
                          << std::to_string(((double) std::chrono::duration_cast<std::chrono::microseconds>(
                                  std::chrono::high_resolution_clock::now() - start).count() / 1000000.0)) << "s"
                          << std::endl << std::endl;
            }
            if (parameters.save) {
                parameters.detailed_evaluation.headerValueInsert(
                        {"Graph", "Size", "Edges", "Parameters", "GeneratorSize", "CoreIterations", "Seed", "Results",
                         "Iteration", "ClosureSize", "CoreSize", "Total Runtime"},
                        {_graph.GetName(), std::to_string(_graph.nodes()), std::to_string(_graph.edges()), "",
                         std::to_string(parameters.generator_size), std::to_string(parameters.core_iterations),
                         std::to_string(i + c_iterations * parameters.seed), "", std::to_string(i),
                         std::to_string(closureParameters.closed_set.size()), std::to_string(overlap.size()),
                         std::to_string(((double) std::chrono::duration_cast<std::chrono::microseconds>(
                                 std::chrono::high_resolution_clock::now() - start).count() / 1000000.0))});
            }
            if (i != 0 && parameters.core_iterations == -1 && parameters.intersection_loss.back() == 0) {
                break;
            }
            ++iteration_count;
        }
        parameters.core_nodes.clear();
        for (auto elem: overlap) {
            parameters.core_nodes.emplace_back(elem);
        }
        //Measure runtime before getting statistics
        parameters.runtime = ((double) std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - start).count() /
                                    1000000.0);
        if (parameters.save) {
            parameters.core_evaluation.headerValueInsert(
                    {"Method", "Graph", "Size", "Edges", "Parameters", "NumIterations", "GrowSteps", "CorePercentage",
                     "GeneratorSize", "Seed", "Results", "CoreSize", "Iterations", "Runtime"},
                    {"RandomCore", _graph.GetName(), std::to_string(_graph.nodes()), std::to_string(_graph.edges()), "",
                     std::to_string(parameters.core_iterations), "", "", std::to_string(parameters.generator_size),
                     std::to_string(parameters.seed), "", std::to_string(parameters.core_nodes.size()),
                     std::to_string(iteration_count), std::to_string(parameters.runtime)});
        }
    }

    GraphStruct& _graph;
};


#endif //TESTGRAPHLIB_COREALGORITHM_H
