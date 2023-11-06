//
// Created by florian on 15.08.23.
//

#ifndef GCOREAPPROXIMATION_COREGROWALGORITHM_H
#define GCOREAPPROXIMATION_COREGROWALGORITHM_H


/**
 * Input parameters for the core grow algorithm introduced by Marc and Subelj in their paper "Convexity in Complex Networks" (2016)
 * The algorithm is a probabilistic algorithm that finds the geodesic core of a graph, in each run it chooses a random start vertex and then grows the geodesic closed set by adding randomly adding a neighbor of the current set and calculating the closure of the new set.
 * num_runs: number of runs of the algorithm
 * grow_steps: number of growth steps per run, i.e. the number of elements added in the iterative growth process
 * core_percentage: percentage of runs in which a node has to be in the core to be considered a core node
 * _seed: _seed for the random number generator
 * _print: if true prints the growth steps and the closure size
 * save: if true saves the core nodes vector
 * output_path: path to the output file
 * core_nodes: vector of core nodes
 * runtime: runtime of the algorithm
 */
struct CoreGrowAlgorithmParameters {
    // Input parameters
    int num_runs = 100;
    int grow_steps = 15;
    double core_percentage = 0.9;
    int seed = 0;
    bool print = false;
    bool save = true;

    // Output parameters
    FileEvaluation core_evaluation;
    FileEvaluation detailed_evaluation;
    std::vector<NodeId> core_nodes = std::vector<NodeId>();
    double runtime = 0.0;
};

class CoreGrowAlgorithm {
public:
    /**
     * Constructor
     * @param graph the graph on which the algorithm is run
     */
    explicit CoreGrowAlgorithm(GraphStruct &graph) : _graph(graph) {}

    /**
     * Runs the core grow algorithm
     * @param parameters
     */
    void Run(CoreGrowAlgorithmParameters& parameters);

private:
    GraphStruct& _graph;
};

void CoreGrowAlgorithm::Run(CoreGrowAlgorithmParameters& parameters){
    // create the _graph closure
    GraphClosure gc = GraphClosure(_graph);
    GraphClosureParameters closureParameters;
    // save the results
    parameters.core_evaluation = FileEvaluation();
    parameters.detailed_evaluation = FileEvaluation();
    // set start time
    auto start = std::chrono::high_resolution_clock::now();
    // set the generator using the _seed
    std::mt19937 generator(parameters.seed);
    // neighbors of the core nodes (in fact the out_edges thus duplicates are possible)
    std::vector<NodeId> core_neighbors = std::vector<NodeId>();
    // vector measuring the size evolution of the core
    std::vector<std::vector<int>> core_size_evolution = std::vector<std::vector<int>>(parameters.num_runs, std::vector<int>(parameters.grow_steps + 1, 0));
    for (auto & x : core_size_evolution) {
        if (!x.empty()){
            x[0] = 1;
        }
    }
    std::vector<NodeId> core_nodes = std::vector<NodeId>(_graph.nodes(), 0);

    // iterate over the number of runs
    for (int i = 0; i < parameters.num_runs; ++i)
    {
        // get the start vertex
        NodeId start_vertex = std::uniform_int_distribution<NodeId>(0, _graph.nodes() - 1)(generator);
        core_neighbors.clear();
        for(auto node : _graph.get_neighbors(start_vertex)){
            core_neighbors.emplace_back(node);
        }

        // get the closure of the start vertex
        std::set<NodeId> start_set = std::set<NodeId>({start_vertex});
        closureParameters.closed_set = start_set;


        // iterate over the growth steps
        for (int j = 0; j < parameters.grow_steps; ++j)
        {
            // _print the growth step
            if (parameters.print)
            {
                std::cout << "Run: " << i << " Grow step: " << j << std::endl;
            }

            // delete all core neighbors that are already in the closureParameters.closed_set
            std::vector<NodeId> core_neighbors_copy = std::vector<NodeId>();
            for (auto node : core_neighbors) {
                if (closureParameters.closed_set.find(node) == closureParameters.closed_set.end())
                {
                    core_neighbors_copy.emplace_back(node);
                }

            }
            core_neighbors.clear();
            for (auto node : core_neighbors_copy) {
                core_neighbors.emplace_back(node);
            }

            // add the neighbors of all added elements to the core neighbors if the neighbor is not already in the core nodes
            for (auto added_element : closureParameters.added_elements)
            {
                // get the neighbors of the added element
                // iterate over the neighbors of the added element
                for (auto added_element_neighbor : _graph.get_neighbors(added_element))
                {
                    // add the neighbor to the core neighbors if it is not already in the core nodes
                    if (closureParameters.closed_set.find(added_element_neighbor) == closureParameters.closed_set.end())
                    {
                        core_neighbors.emplace_back(added_element_neighbor);
                    }
                }
            }
            if (!core_neighbors.empty()) {
                // get a random element from the core neighbors
                NodeId random_element = core_neighbors[std::uniform_int_distribution<NodeId>(0,
                                                                                             (NodeId) core_neighbors.size() -
                                                                                             1)(generator)];
                // calculate the closure of the core + the random element
                closureParameters.input_set = closureParameters.closed_set;
                closureParameters.element_to_add = random_element;
                gc.closure(closureParameters);
                // add the random element to the added elements
                closureParameters.added_elements.insert(random_element);

                core_size_evolution[i][j + 1] = (int) closureParameters.added_elements.size();
            }
        }
        // _print the closure size
        if (parameters.print)
        {
            std::cout << "Closure size: " << closureParameters.closed_set.size() << std::endl;
        }
        // add the core nodes to the output parameters
        for (auto node : closureParameters.closed_set)
        {
            core_nodes[node] += 1;
        }
        // save the results
        parameters.detailed_evaluation.headerValueInsert({"Graph", "Size", "Edges", "Parameters", "NumRuns", "GrowSteps", "CorePercentage", "Seed", "Results", "Iteration", "ClosureSize", "CoreSize", "Total Runtime"},
                                             {_graph.GetName(), std::to_string(_graph.nodes()), std::to_string(_graph.edges()), "", std::to_string(parameters.num_runs), std::to_string(parameters.grow_steps), std::to_string(parameters.core_percentage), std::to_string(parameters.seed), "", std::to_string(i), std::to_string(closureParameters.closed_set.size()), "", std::to_string(parameters.runtime)});
    }
    // _print the core nodes vector
    if (parameters.print)
    {
        std::cout << "Core nodes: " << std::endl;
        for (int i = 0; i < core_nodes.size(); ++i)
        {
            std::cout << i << ": " << core_nodes[i] << std::endl;
        }
    }
    // iterate over core_nodes and add to output parameters if the node is in core node percentage of the runs
    for (int i = 0; i < core_nodes.size(); ++i)
    {
        if (core_nodes[i] >= int(parameters.core_percentage * parameters.num_runs + 0.99))
        {
            parameters.core_nodes.emplace_back(i);
        }
    }
    // set the end time
    parameters.runtime = (double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1e6;

    parameters.detailed_evaluation.headerValueInsert({"Graph", "Size", "Edges", "Parameters", "NumRuns", "GrowSteps", "CorePercentage", "Seed", "Results", "Iteration", "ClosureSize", "CoreSize", "Runtime"},
                                         {_graph.GetName(), std::to_string(_graph.nodes()), std::to_string(_graph.edges()), "", std::to_string(parameters.num_runs), std::to_string(parameters.grow_steps), std::to_string(parameters.core_percentage), std::to_string(parameters.seed), "", "", "", std::to_string(parameters.core_nodes.size()), std::to_string(parameters.runtime)});


    parameters.core_evaluation.headerValueInsert({"Method", "Graph", "Size", "Edges", "Parameters", "NumIterations", "GrowSteps", "CorePercentage", "GeneratorSize", "Seed", "Results", "CoreSize", "Iterations", "Runtime"},
                                     {"GrowCore", _graph.GetName(), std::to_string(_graph.nodes()), std::to_string(_graph.edges()), "", std::to_string(parameters.num_runs), std::to_string(parameters.grow_steps), std::to_string(parameters.core_percentage), "", std::to_string(parameters.seed), "", std::to_string(parameters.core_nodes.size()), "", std::to_string(parameters.runtime)});
}

#endif //GCOREAPPROXIMATION_COREGROWALGORITHM_H
