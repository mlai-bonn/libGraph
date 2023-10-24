//
// Created by Florian on 09.03.2021.
//


#ifndef HOPS_HOPS_H
#define HOPS_HOPS_H

#include <utility>
#include <vector>
#include <set>
#include <string>
#include "Evaluation.h"
#include "RootedPattern.h"
#include "HopsParameters.h"
#include <omp.h>
#include <list>
#include <unordered_set>
#include "../io/FileEvaluation.h"
#include "GraphDataStructures/GraphBase.h"
#include "GraphDataStructures/GraphDirectedBase.h"


class Hops {
public:
    /// Default constructor
    Hops() = default;
    /// Constructor which gets the graph data and some parameters for running the approximation
    /// @param graphs input graphs
    /// @param parameters input parameters such as Save path etc.
    explicit Hops(GraphData<GraphStruct>& graphs, const std::string& out_path = "") : graphs(graphs), out_path(out_path) {
        for (auto & graph : this->graphs.graphData) {
            Nodes order;
            //GraphStruct tree;
            //GraphStruct::DFS(graph, tree, order);
            //GraphStruct::ReorderGraph(graph, order);
        }
    };
    /// Constructor which gets the graph data and some parameters for running the approximation
    /// @param graphs input graphs
    /// @param parameters input parameters such as Save path etc.
    explicit Hops(const std::string& in_path, const std::string& out_path = "") : out_path(out_path) {
        this->graphs = GraphData<GraphStruct>(in_path);
        for (auto & graph : this->graphs.graphData) {
            Nodes order;
            //GraphStruct tree;
            //GraphStruct::DFS(graph, tree, order);
            //GraphStruct::ReorderGraph(graph, order);
        }
    };

    /// Main function of the hops class runs the hops algorithm with different run parameters
    /// @param graphId id of the graph used for estimation
    /// @param pattern Pattern used for estimation
    /// @param runParameters Parameters for the run of hops e.g. number of iterations, fixed runtime, run labeled or not etc...
    void Run(size_t graphId, GraphStruct& pattern, RunParameters rParameters = RunParameters());
    /// Estimate automorphisms for a certain graph
    /// \param graphId
    /// \param rParameters
    void Automorphisms(size_t graphId, RunParameters rParameters = RunParameters());

    const GraphData<GraphStruct>& GetGraphs(){return graphs;};

    struct RunProps;

private:
    /// Preprocess the pattern, mainly create rooted pattern from some pattern spanning tree and returns the label type for the run
    /// @param spanTree id of the graph used for estimation
    /// @param pattern_seed seed for some randomness of the spanning tree of the pattern
    LABEL_TYPE PatternPreprocessing(const DGraphStruct & spanTree, int pattern_seed = 0);
    void InitEstimation(unsigned int& id, RunProps& runProps, const std::vector<NodeId>& PatternRootImages, bool random = true);
    /// Main function for the estimation of how often an unlabeled pattern occurs in the big graph
    /// @param id id of the estimation
    /// @param accumulatedMean accumulated estimation number for all estimations up to now
    /// @param estimationVector vector of estimations
    /// @param zeroIterations counter for the iterations which are zero
    /// @param runProps extended properties used in the run
    bool UnlabeledEmbedding(unsigned int id, long double& accumulatedMean, long double& accumulatedDeviation, std::vector<UInt64>& estimationVector, UInt64& zeroIterations, RunProps& runProps);
    /// Main function for the estimation of how often a labelled pattern occurs in the big graph
    /// @param id id of the estimation
    /// @param accumulatedEstimation accumulated estimation number for all estimations up to now
    /// @param estimationVector vector of estimations
    /// @param zeroIterations counter for the iterations which are zero
    /// @param runProps extended properties used in the run
    bool LabeledEmbeddings(unsigned int id, long double& accumulatedEstimation, long double& accumulatedDeviation, std::vector<UInt64>& estimationVector, UInt64& zeroIterations, RunProps& runProps);
    bool UnlabeledGraphEditDistance(unsigned int id, unsigned int &currentEditDistance, RunProps &runProps);
    bool LabeledGraphEditDistance(unsigned int id, unsigned int& currentEditDistance, RunProps &runProps);
    bool PickRandomNeighbors(unsigned int id, unsigned int sourceSize, unsigned int targetSize, NodeId currentNodeId, const Nodes& currentNodeNeighbors, unsigned int& possibleNeighbors, const Nodes &orderNodes, RunProps& runProps, int label = -1) const;
    bool FindEmbedding(unsigned int id, unsigned int sourceSize, unsigned int targetSize, NodeId currentNodeId, unsigned int& numberOfEmbeddedNodes, unsigned int &numberOfEmbeddedEdges, const Nodes &orderNodes, RunProps& runProps, int label = -1) const;
    bool CheckValidGraphEmbedding(const Nodes& nodes, const RunProps& runProps) const;
    bool CheckValidGraphEmbedding(const Nodes& nodes, unsigned int& embeddedEdges, const RunProps& runProps) const;
    NodeId getRootNodeByCondition(ROOT_NODE_CONDITION condition, std::mt19937_64& gen);
    void EvaluateResult(std::vector<UInt64>& estimations, std::vector<long double>& snapshots);
    void EvaluateResult(std::vector<long double>& snapshots);
    void EvaluateResult(int approximatedGED, HOPS_TYPE hopsType = HOPS_TYPE::ESTIMATION);


    GraphData<GraphStruct> graphs;
    GraphStruct* currentGraph = nullptr;
    GraphStruct* currentPattern = nullptr;

    //If CurrentPattern is not a tree this vector contains all edges not in
    std::vector<EDGES> nonTreeEdges;
    DGraphStruct spanningTree;
    INDEX spanningTreeRootNode = std::numeric_limits<INDEX>::max();
    RootedPattern rootedPattern;
    Nodes possibleGraphImagesOfPatternRoot;
    std::uniform_int_distribution<NodeId> initialImageDistribution;
    UInt64 logValueMultiplier = 100000000;
    int currentGraphMaxDegree = 0;
    bool currentPatternIsTree = false;
    std::vector<UInt64> preComputedLogValues;
    std::string out_path;
    RunParameters runParameters;
public:

    Evaluation hopsEvaluation = Evaluation();

    struct RunProps{
        RunProps(const GraphStruct* graph, const GraphStruct* pattern, int graphMaxDegree, LABEL_TYPE labelType = LABEL_TYPE::UNLABELED, RUN_TYPE type = RUN_TYPE::VECTOR, int seed = 0) :
                treeGraphMap(Nodes(pattern->nodes(), std::numeric_limits<INDEX>::max())),
                NewImagePositions(std::vector<NodeId>(pattern->nodes(), 0)),
                runType(type),
                labelType(labelType),
                randomNeighbors(std::vector<NodeId>(graphMaxDegree, 0)),
                randomNeighborsSwapPairs(std::vector<std::pair<NodeId, NodeId>>(graphMaxDegree, {0,0})),
                graphImages(std::unordered_set<NodeId>()),
                savedNeighbors(Nodes()),
                gen(std::mt19937_64(seed)),
                seed(seed){
            std::iota(randomNeighbors.begin(), randomNeighbors.end(), 0);
            if (type == RUN_TYPE::VECTOR) {
                savedNeighbors = Nodes(graph->nodes(), 0);
            }
        };

        bool VisitNeighbor(unsigned int Id, NodeId neighborNodeId);

        INDEX currentRootNode = 0;
        INDEX currentNodeIteration = 0;
        INDEX counter = 0;
        int current_iteration_step = 0;
        std::vector<long double> estimation_snapshots;
        std::vector<UInt64> iteration_snapshots;

        Nodes treeGraphMap;
        std::vector<NodeId> NewImagePositions;
        std::vector<NodeId> randomNeighbors;
        std::vector<std::pair<NodeId, NodeId>> randomNeighborsSwapPairs;
        unsigned int swapCount = 0;
        Nodes savedNeighbors;
        std::unordered_set<NodeId> graphImages;
        std::mt19937_64 gen;
        int seed;
        RUN_TYPE runType;
        LABEL_TYPE labelType;

        void swapNeighbors(unsigned int idx_a, unsigned int idx_b);

        void swapBackNeighbors();

        bool IsNeighborVisited(unsigned int Id, NodeId neighborNodeId);

        std::vector<UInt64> privateEstimations;
        std::vector<long double> privateSnapshots;

        std::chrono::time_point<std::chrono::high_resolution_clock , std::chrono::duration<double>> currentTime;
        std::chrono::time_point<std::chrono::high_resolution_clock, std::chrono::duration<double>> endTime;
        std::chrono::time_point<std::chrono::high_resolution_clock, std::chrono::duration<double>> lastSnapShotTime;
        bool runAlgorithm = false;



    };

    //void CheckEstimationFinished(RunProps &props, UInt64 OverallIterations, const Nodes& possibleRootImages);

    void MergeEstimationValues(RunProps& runProps, std::vector<UInt64>& estimations, std::vector<long double>& snapShotEstimation) const;


    void FastRun(long double accumulatedMean, long double accumulatedStd,
                 long double nodeMean, long double nodeStd, UInt64 OverallIterations, UInt64 OverallZeroIterations,
                 unsigned int accumulatedApproximatedGED, unsigned int approximatedGED,
                 std::vector<UInt64> &estimations,
                 std::vector<long double> &snapShotEstimation, RunProps &runProps, unsigned int Id);
    void RandomRun(long double accumulatedMean, long double accumulatedStd,
                   long double nodeMean, long double nodeStd, UInt64 OverallIterations, UInt64 OverallZeroIterations,
                   unsigned int accumulatedApproximatedGED, unsigned int approximatedGED,
                   std::vector<UInt64> &estimations,
                   std::vector<long double> &snapShotEstimation, RunProps &runProps, unsigned int Id);
};


inline NodeId Hops::getRootNodeByCondition(ROOT_NODE_CONDITION condition, std::mt19937_64& gen){
    NodeId rootNodeId = 0;
    size_t MaxDegree = 0;
    size_t MinDegree = std::numeric_limits<size_t>::max();
    switch (condition) {
        case TREE_GIVEN:
            for (NodeId nodeId = 0; nodeId < this->spanningTree.nodes(); ++nodeId) {
                if (this->spanningTree.degree(nodeId) == 0){
                    return nodeId;
                }
            }
            break;
        case NODE_GIVEN:
            return this->spanningTreeRootNode;
            break;
        case MAX_DEGREE:
            //Get Node with max degree as root node
            for (NodeId nodeId = 0; nodeId < this->currentPattern->nodes(); ++nodeId){
                INDEX Degree = this->currentPattern->degree(nodeId);
                if (Degree > MaxDegree) {
                    MaxDegree = Degree;
                    rootNodeId = nodeId;
                }
            }
            return rootNodeId;
            break;
        case MIN_DEGREE:
            //Get Node with max degree as root node
            for (NodeId nodeId = 0; nodeId < this->currentPattern->nodes(); ++nodeId){
                INDEX Degree = this->currentPattern->degree(nodeId);
                if (Degree < MinDegree) {
                    MinDegree = Degree;
                    rootNodeId = nodeId;
                }
            }
            return rootNodeId;
            break;
        case MIN_LABEL_BIG_GRAPH:
            for(const std::pair<Label, int> LabelFrequency : this->currentGraph->labelFrequencyMap){
                if (this->currentPattern->labelMap.find(LabelFrequency.first) != this->currentPattern->labelMap.end()){
                    int RandNumber = std::uniform_int_distribution<>(0, (int) this->currentPattern->labelMap[LabelFrequency.first].size() - 1)(gen);
                    return this->currentPattern->labelMap[LabelFrequency.first][RandNumber];
                }
            }
            break;
    }
    return 0;
}

inline LABEL_TYPE Hops::PatternPreprocessing(const DGraphStruct & spanTree, int pattern_seed) {
    this->currentPatternIsTree = this->currentPattern->IsTree();

    LABEL_TYPE labelType = this->currentGraph->labelType;
    if (this->runParameters.labelType == LABEL_TYPE::UNLABELED || this->currentGraph->labelType == LABEL_TYPE::UNLABELED || this->runParameters.labelType != this->currentGraph->labelType)
    {
        labelType = LABEL_TYPE::UNLABELED;
    }

    NodeId PatternRootNode = std::numeric_limits<NodeId>::max();
    std::mt19937_64 randomGenerator(pattern_seed);
    //Generate or check spanning tree for patternId
    if (spanTree.nodes() > 0 && !this->currentPattern->CheckSpanningTree(spanTree)){
        throw std::invalid_argument("Spanning tree is not valid");
    }
    else{
        //GetRootNode of tree
        //Spanning tree is already given
        if (spanTree.nodes() > 0){
            this->spanningTree = DGraphStruct(spanTree);

            PatternRootNode = getRootNodeByCondition(TREE_GIVEN, randomGenerator);
        }
            //Spanning tree is not given
        else {
            //Root node is given
            if (this->spanningTreeRootNode != std::numeric_limits<INDEX>::max()) {
                PatternRootNode = getRootNodeByCondition(NODE_GIVEN, randomGenerator);
                this->spanningTree = DGraphStruct::GetBFSTree(*this->currentPattern, PatternRootNode);
            }
                // Get some root node by some condition
            else{
                ROOT_NODE_CONDITION Condition;
                if (labelType == LABEL_TYPE::UNLABELED){
                    Condition = MAX_DEGREE;
                }
                else{
                    Condition = MIN_LABEL_BIG_GRAPH;
                }
                PatternRootNode = getRootNodeByCondition(Condition, randomGenerator);
                this->spanningTree = DGraphStruct::GetBFSTree(*this->currentPattern, PatternRootNode);
            }

        }
        //Generate the rooted pattern from the pattern graph
        this->rootedPattern = RootedPattern(this->spanningTree, PatternRootNode, this->runParameters.labelType);

        if (!this->currentPattern->isTree) {
            this->nonTreeEdges = std::vector<EDGES>(this->currentPattern->nodes(), EDGES());
            for (NodeId nodeId = 0; nodeId < this->currentPattern->nodes(); ++nodeId) {
                NodeId Source = nodeId;
                for (int neighborIndex = 0; neighborIndex < this->currentPattern->degree(nodeId); ++neighborIndex) {
                    NodeId Destination = this->currentPattern->neighbor(nodeId, neighborIndex);
                    NodeId RootedPatternSource = this->rootedPattern.GetBFSOrderIndexByNodeId(Source);
                    NodeId RootedPatternDestination = this->rootedPattern.GetBFSOrderIndexByNodeId(Destination);
                    if (!this->spanningTree.edge(Source, Destination, false)) {
                        this->nonTreeEdges[std::max(RootedPatternSource, RootedPatternDestination)].emplace_back(RootedPatternSource, RootedPatternDestination);
                    }
                }
            }
        }
        //Define possible root node mappings (using the degree of the root node)
        this->possibleGraphImagesOfPatternRoot.clear();
        if (labelType == LABEL_TYPE::UNLABELED){
            for (NodeId i = 0; i < this->currentGraph->nodes(); ++i) {
                if(this->currentGraph->degree(i) >= this->currentPattern->degree(PatternRootNode)){
                    this->possibleGraphImagesOfPatternRoot.emplace_back(i);
                }
            }
        }
        else{
            Label RootNodeLabel = this->currentPattern->label(PatternRootNode);
            this->possibleGraphImagesOfPatternRoot = this->currentGraph->labelMap[RootNodeLabel];
        }
    }

    return labelType;
}



inline void Hops::Run(size_t graphId, GraphStruct& pattern, RunParameters rParameters) {
    this->runParameters = rParameters;
    if (this->runParameters.thread_num == -1){
        this->runParameters.thread_num = omp_get_max_threads();
    }
    omp_set_num_threads(this->runParameters.thread_num);
    auto start = std::chrono::high_resolution_clock::now();
    this->currentPattern = &pattern;

    //Value of accumulated estimation and counters for overall (zero) iterations
    long double accumulatedMean = 0;
    long double accumulatedStd = 0;

    long double nodeMean = 0;
    long double nodeStd = 0;
    //Counter of overall iterations
    UInt64 OverallIterations = 0;
    //Counter for overall zero iterations
    UInt64 OverallZeroIterations = 0;

    //For GED approximation
    unsigned int accumulatedApproximatedGED = std::numeric_limits<unsigned int>::max();
    unsigned int approximatedGED = std::numeric_limits<unsigned int>::max();

    //Vectors for estimations values
    std::vector<UInt64> estimations;
    std::vector<long double> snapShotEstimation;

    //Preprocess the graph information, i.e. precompute the possible combination values by considering the max degree in the graph
    if (this->currentGraph != &graphs[graphId]) {
        this->currentGraph = &graphs[graphId];
        this->currentGraphMaxDegree = this->currentGraph->maxDegree;
        //Recalculate the precomputed log values
        this->preComputedLogValues = std::vector<UInt64>(this->currentGraphMaxDegree + 1, 0);
        long double Value = 0;
        this->preComputedLogValues[0] = std::numeric_limits<UInt64>::max();
        for (int i = 2; i <= this->currentGraphMaxDegree; ++i) {
            Value += log(i);
            if (Value > (double) std::numeric_limits<UInt64>::max()/ (double) this->logValueMultiplier){
                throw std::overflow_error("Max Degree is too big!!!");
            }
            this->preComputedLogValues[i] = (UInt64) (Value * this->logValueMultiplier);
        }

    }
    this->spanningTreeRootNode = this->runParameters.spanningTreeRootNode;
    //Preprocess the pattern graph
    LABEL_TYPE labelType = this->PatternPreprocessing(this->runParameters.spanningTree, this->runParameters.seed);

    //Set the run properties
    RunProps runProps = RunProps(this->currentGraph, this->currentPattern, this->currentGraphMaxDegree, labelType, this->runParameters.run_type, this->runParameters.seed);

    //
    if (!this->runParameters.runtime.empty()){
        this->runParameters.iteration_per_node = {0};
    }
    this->hopsEvaluation = Evaluation(*this->currentGraph, *this->currentPattern, this->runParameters);
    //save hops preprocessing runtime
    this->hopsEvaluation.preprocessingTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);



    //runProps.endTime = start + std::chrono::microseconds ((std::size_t) (this->runParameters.runtime*1e6));

    //set run id, this will be used to determine the visited nodes by comparing ids
    unsigned int Id = 0;


    if(!this->runParameters.runtime.empty())
    {
        this->runParameters.iteration_per_node = {0};
        snapShotEstimation.resize(this->runParameters.runtime.size());
        //Run the estimation for random nodes in the graph
        RandomRun(accumulatedMean, accumulatedStd, nodeMean, nodeStd, OverallIterations, OverallZeroIterations,
                  accumulatedApproximatedGED, approximatedGED, estimations, snapShotEstimation, runProps, Id);
    }
    else {
        snapShotEstimation.resize(this->runParameters.iteration_per_node.size());
        //Run the estimation for each node in the graph
        FastRun(accumulatedMean, accumulatedStd, nodeMean, nodeStd, OverallIterations, OverallZeroIterations,
                accumulatedApproximatedGED, approximatedGED, estimations, snapShotEstimation, runProps, Id);
    }
}

void Hops::FastRun(long double accumulatedMean, long double accumulatedStd,
                   long double nodeMean, long double nodeStd, UInt64 OverallIterations, UInt64 OverallZeroIterations,
                   unsigned int accumulatedApproximatedGED, unsigned int approximatedGED,
                   std::vector<UInt64> &estimations, std::vector<long double> &snapShotEstimation,
                   Hops::RunProps &runProps, unsigned int Id) {
    //start clock for hops runtime measure
    auto start = std::chrono::high_resolution_clock::now();
    ////run the estimation in parallel mode
#pragma omp parallel default(none) firstprivate(runProps, approximatedGED, nodeMean, nodeStd, Id) shared(estimations, snapShotEstimation, accumulatedApproximatedGED, possibleGraphImagesOfPatternRoot) reduction(+: accumulatedMean, accumulatedStd, OverallIterations, OverallZeroIterations)
    {
        int thread_num = omp_get_thread_num();
        runProps.runAlgorithm = true;
        runProps.currentTime = std::chrono::high_resolution_clock::now();
        runProps.lastSnapShotTime = std::chrono::high_resolution_clock::now() - std::chrono::hours(1);
        runProps.estimation_snapshots.resize(runParameters.iteration_per_node.size());

#pragma omp for schedule(dynamic)
        for (auto currentNode : possibleGraphImagesOfPatternRoot) {
            //set current node
            runProps.currentRootNode = currentNode;
            nodeMean = 0;
            nodeStd = 0;
            //reset iteration to zero for the
            runProps.current_iteration_step = 0;
            //set the seed for this node
            runProps.gen.seed(runProps.seed + currentNode);
            for (int iter = 0; iter < runParameters.iteration_per_node.back(); ++iter) {
                //get estimation of one embedding
                InitEstimation(Id, runProps, possibleGraphImagesOfPatternRoot, false);
                switch (runParameters.hops_type) {
                    //Hops embedding algorithm
                    case HOPS_TYPE::ESTIMATION:
                        //Start the embedding
                        //Unlabeled Case
                        if (runProps.labelType == LABEL_TYPE::UNLABELED) {
                            UnlabeledEmbedding(Id, nodeMean, nodeStd, runProps.privateEstimations,
                                               OverallZeroIterations, runProps);
                        }
                            //Labeled Case
                        else {
                            LabeledEmbeddings(Id, nodeMean, nodeStd, runProps.privateEstimations, OverallZeroIterations,
                                              runProps);
                        }
                        ++OverallIterations;
                        break;
                        //TODO add graph edit distance algorithm
                    case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
                        unsigned int approxGED = static_cast<int>(currentPattern->nodes() + currentPattern->edges());
                        runProps.labelType = LABEL_TYPE::UNLABELED;
                        UnlabeledGraphEditDistance(Id, approxGED, runProps);
                        if (runProps.labelType != LABEL_TYPE::UNLABELED) {
                            for (int i = 0; i < runProps.treeGraphMap.size(); ++i) {
                                int TreeNode = i;
                                NodeId GraphNode = runProps.treeGraphMap[i];
                                if(GraphNode >= 0){
                                    if(currentPattern->label(TreeNode) != currentGraph->label(GraphNode)){
                                        ++approxGED;
                                    }
                                }
                            }
                        }
                        ++OverallIterations;
                        approximatedGED = std::min(approximatedGED, approxGED);
                        if (approxGED == 0){
                            runProps.runAlgorithm = false;
                        }
                        else {
                            //CheckEstimationFinished(runProps, OverallIterations, possibleGraphImagesOfPatternRoot);
                        }
                        break;
                }
                //Collect on the fly results
                if (iter == runParameters.iteration_per_node[runProps.current_iteration_step] - 1){
                    runProps.estimation_snapshots[runProps.current_iteration_step] += nodeMean;
                    ++runProps.current_iteration_step;
                }
            }
            accumulatedMean += nodeMean;
            accumulatedStd += nodeStd;
        }

        if (runParameters.iteration_per_node.size() > 1){
            for (int i = 0; i < runProps.estimation_snapshots.size(); ++i) {
                runProps.estimation_snapshots[i] /= (possibleGraphImagesOfPatternRoot.size() * runParameters.iteration_per_node[i]);
            }
        }

        if (runParameters.snapshot_time != -1) {
            runProps.privateSnapshots.emplace_back((double) accumulatedMean / (double) OverallIterations);
        }

        //merge evaluation values from different threads in one large vector if needed (do this only with one thread)
        if (runParameters.iteration_per_node.size() > 1) {
#pragma omp critical
            {
                for (int i = 0; i < snapShotEstimation.size(); ++i) {
                    snapShotEstimation[i] += runProps.estimation_snapshots[i];
                }
            }
        }

        if (!runParameters.single_number || runParameters.snapshot_time != -1) {
#pragma omp critical
            {
                switch (runParameters.hops_type) {
                    case HOPS_TYPE::ESTIMATION:
                        MergeEstimationValues(runProps, estimations, snapShotEstimation);
                        break;
                    case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
                        accumulatedApproximatedGED = std::min(accumulatedApproximatedGED, approximatedGED);
                        break;
                }
            }
        }
    }

    //measure hops time
    hopsEvaluation.hopsRuntime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
    hopsEvaluation.hopsIterations = OverallIterations;

    //evaluate the estimation result of hops
    switch (runParameters.hops_type) {
        case HOPS_TYPE::ESTIMATION:
            //save hops runtime and estimation values and evaluate
            if (runParameters.single_number){
                hopsEvaluation.hopsEstimation = accumulatedMean;
                hopsEvaluation.hopsStd = accumulatedStd;
                hopsEvaluation.hopsZeroIterations = OverallZeroIterations;
                EvaluateResult(snapShotEstimation);
            }else {
                hopsEvaluation.hopsStd = accumulatedStd;
                hopsEvaluation.hopsZeroIterations = OverallZeroIterations;
                EvaluateResult(estimations, snapShotEstimation);
            }
            break;
        case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
            EvaluateResult(accumulatedApproximatedGED, HOPS_TYPE::GRAPH_EDIT_DISTANCE);
            break;
    }
}

void Hops::RandomRun(long double accumulatedMean, long double accumulatedStd,
                   long double nodeMean, long double nodeStd, UInt64 OverallIterations, UInt64 OverallZeroIterations,
                   unsigned int accumulatedApproximatedGED, unsigned int approximatedGED,
                   std::vector<UInt64> &estimations, std::vector<long double> &snapShotEstimation,
                   Hops::RunProps &runProps, unsigned int Id) {
    //create vector for snapshot times
    std::vector<std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds>> snapshot_times;
    //start clock for hops runtime measure
    auto start = std::chrono::high_resolution_clock::now();
    runProps.endTime = start + std::chrono::microseconds ((std::size_t) (this->runParameters.runtime.back()*1e6));
    // determine snapshot times using runParameters.runtime and start
    for (double i : runParameters.runtime) {
        snapshot_times.emplace_back(start + std::chrono::microseconds ((std::size_t) (i*1e6)));
    }

    std::vector<UInt64> snapShotIterations;
    snapShotIterations.resize(runParameters.runtime.size());

    ////run the estimation in parallel mode
#pragma omp parallel default(none) firstprivate(start, runProps, approximatedGED, nodeMean, nodeStd, Id) shared(snapshot_times, estimations, snapShotEstimation, snapShotIterations, accumulatedApproximatedGED, possibleGraphImagesOfPatternRoot) reduction(+: accumulatedMean, accumulatedStd, OverallIterations, OverallZeroIterations)
    {
        int thread_num = omp_get_thread_num();
        runProps.runAlgorithm = true;
        runProps.currentTime = std::chrono::high_resolution_clock::now();
        runProps.lastSnapShotTime = std::chrono::high_resolution_clock::now() - std::chrono::hours(1);
        runProps.estimation_snapshots.resize(runParameters.runtime.size());
        runProps.current_iteration_step = 0;



        // while time is not greater than runProps.endTime
        runProps.currentTime = std::chrono::high_resolution_clock::now();

        // get thread id
        auto threadId = omp_get_thread_num();
        runProps.gen.seed(runProps.seed + threadId);
        while (runProps.current_iteration_step < runParameters.runtime.size()) {
            //get estimation of one embedding
            InitEstimation(Id, runProps, possibleGraphImagesOfPatternRoot, true);
            switch (runParameters.hops_type) {
                //Hops embedding algorithm
                case HOPS_TYPE::ESTIMATION:
                    //Start the embedding
                    //Unlabeled Case
                    if (runProps.labelType == LABEL_TYPE::UNLABELED) {
                        UnlabeledEmbedding(Id, accumulatedMean, accumulatedStd, runProps.privateEstimations,
                                           OverallZeroIterations, runProps);
                    }
                        //Labeled Case
                    else {
                        LabeledEmbeddings(Id, accumulatedMean, accumulatedStd, runProps.privateEstimations,
                                          OverallZeroIterations,
                                          runProps);
                    }
                    ++OverallIterations;
                    break;
                    //TODO add graph edit distance algorithm
                case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
                    unsigned int approxGED = static_cast<int>(currentPattern->nodes() + currentPattern->edges());
                    runProps.labelType = LABEL_TYPE::UNLABELED;
                    UnlabeledGraphEditDistance(Id, approxGED, runProps);
                    if (runProps.labelType != LABEL_TYPE::UNLABELED) {
                        for (int i = 0; i < runProps.treeGraphMap.size(); ++i) {
                            int TreeNode = i;
                            NodeId GraphNode = runProps.treeGraphMap[i];
                            if (GraphNode >= 0) {
                                if (currentPattern->label(TreeNode) != currentGraph->label(GraphNode)) {
                                    ++approxGED;
                                }
                            }
                        }
                    }
                    ++OverallIterations;
                    approximatedGED = std::min(approximatedGED, approxGED);
                    if (approxGED == 0) {
                        runProps.runAlgorithm = false;
                    } else {
                        //CheckEstimationFinished(runProps, OverallIterations, possibleGraphImagesOfPatternRoot);
                    }
                    break;
            }
            //get current time
            runProps.currentTime = std::chrono::high_resolution_clock::now();
            //Collect on the fly results
            if (runProps.currentTime > snapshot_times[runProps.current_iteration_step]) {
                runProps.estimation_snapshots[runProps.current_iteration_step] = (double) accumulatedMean / (double) OverallIterations;
                runProps.iteration_snapshots.emplace_back(OverallIterations);
                ++runProps.current_iteration_step;
            }
//            //Collect on the fly results
//            if (iter == runParameters.iteration_per_node[runProps.current_iteration_step] - 1){
//                runProps.estimation_snapshots[runProps.current_iteration_step] += nodeMean;
//                ++runProps.current_iteration_step;
//            }
        }

        if (runParameters.snapshot_time != -1) {
            runProps.privateSnapshots.emplace_back((double) accumulatedMean / (double) OverallIterations);
        }

        //merge evaluation values from different threads in one large vector if needed (do this only with one thread)
        if (runParameters.runtime.size() > 1) {
#pragma omp critical
            {
                for (int i = 0; i < snapShotEstimation.size(); ++i) {
                    snapShotEstimation[i] += runProps.estimation_snapshots[i]/(double) this->runParameters.thread_num;
                    snapShotIterations[i] += runProps.iteration_snapshots[i];
                }
            }
        }

        if (!runParameters.single_number || runParameters.snapshot_time != -1) {
#pragma omp critical
            {
                switch (runParameters.hops_type) {
                    case HOPS_TYPE::ESTIMATION:
                        MergeEstimationValues(runProps, estimations, snapShotEstimation);
                        break;
                    case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
                        accumulatedApproximatedGED = std::min(accumulatedApproximatedGED, approximatedGED);
                        break;
                }
            }
        }
    }

    //measure hops time
    hopsEvaluation.hopsRuntime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
    hopsEvaluation.hopsIterations = OverallIterations;
    hopsEvaluation.snapshotIterations = snapShotIterations;

    //evaluate the estimation result of hops
    switch (runParameters.hops_type) {
        case HOPS_TYPE::ESTIMATION:
            //save hops runtime and estimation values and evaluate
            if (runParameters.single_number){
                hopsEvaluation.hopsEstimation = accumulatedMean;
                hopsEvaluation.hopsStd = accumulatedStd;
                hopsEvaluation.hopsZeroIterations = OverallZeroIterations;
                EvaluateResult(snapShotEstimation);
            }else {
                EvaluateResult(estimations, snapShotEstimation);
            }
            break;
        case HOPS_TYPE::GRAPH_EDIT_DISTANCE:
            EvaluateResult(accumulatedApproximatedGED, HOPS_TYPE::GRAPH_EDIT_DISTANCE);
            break;
    }

}

//inline void Hops::CheckEstimationFinished(Hops::RunProps &runProps, UInt64 OverallIterations, const Nodes& possibleRootImages)  {
//    if (this->runParameters.runtime != 0 && runProps.first_run) {
//        if (runProps.first_run && runProps.currentRootNode >= possibleRootImages.size() - 1){
//            auto count =  std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - runProps.first_run_start).count();
//            runProps.first_run_time = (double) count/1e6;
//            std::cout << runProps.first_run_time << std::endl;
//            this->runParameters.iteration_per_node.back() = (int) (this->runParameters.runtime / runProps.first_run_time*1.5);
//            this->hopsEvaluation.parameters.iteration_per_node.back() = this->runParameters.iteration_per_node.back() + 1;
//            std::cout << this->runParameters.iteration_per_node.back() << std::endl;
//            if (this->runParameters.iteration_per_node.back() <= 0){
//                runProps.runAlgorithm = false;
//            }
//            runProps.currentRootNode = 0;
//            runProps.currentNodeIteration = 0;
//            runProps.first_run = false;
//        }
//        else{
//            ++runProps.currentRootNode;
//        }
////        if (OverallIterations % 1000 == 0) {
////            runProps.currentTime = std::chrono::high_resolution_clock::now();
////            if (runProps.endTime < runProps.currentTime) {
////                runProps.runAlgorithm = false;
////            }
////        }
//    }
//    else{
//        INDEX size = possibleRootImages.size();
//        if (runProps.counter == size * this->runParameters.iteration_per_node.back()){
//            runProps.runAlgorithm = false;
//        }
//    }
//}

inline void Hops::MergeEstimationValues(RunProps& runProps, std::vector<UInt64>& estimations, std::vector<long double>& snapShotEstimation) const {
    if (!this->runParameters.single_number) {
        estimations.insert(estimations.end(), runProps.privateEstimations.begin(), runProps.privateEstimations.end());
    }
    if (this->runParameters.snapshot_time != -1) {
        if (snapShotEstimation.empty()) {
            for (long double estimation: runProps.privateSnapshots) {
                snapShotEstimation.emplace_back(estimation / (double) this->runParameters.thread_num);
            }

        } else {
            for (size_t i = 0; i < std::min(snapShotEstimation.size(), runProps.privateSnapshots.size()); ++i) {
                snapShotEstimation[i] += runProps.privateSnapshots[i] / (double) this->runParameters.thread_num;
            }
        }
    }
}

inline void Hops::EvaluateResult(std::vector<long double>& snapshots) {
    auto start = std::chrono::high_resolution_clock::now();
    long double startEstimation = (double) this->possibleGraphImagesOfPatternRoot.size();
    this->hopsEvaluation.hopsEstimation *= startEstimation / (long double) this->hopsEvaluation.hopsIterations;
    this->hopsEvaluation.hopsStd *= (startEstimation * startEstimation) / (long double) this->hopsEvaluation.hopsIterations;
    this->hopsEvaluation.snapshots = snapshots;
    for (long double & snapshot : this->hopsEvaluation.snapshots) {
        snapshot *= startEstimation;
    }

    //Evaluation Runtime
    this->hopsEvaluation.evaluationRuntime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);

    //Print and Save the results
    if (this->runParameters.print){
        this->hopsEvaluation.print();
    }
    if(!this->out_path.empty()){
        this->hopsEvaluation.save(this->out_path);
    }
}

inline void Hops::EvaluateResult(std::vector<UInt64>& estimations, std::vector<long double>& snapshots) {
    auto start = std::chrono::high_resolution_clock::now();
    //Get the results
    this->hopsEvaluation.hopsEstimation = 0;
    this->hopsEvaluation.hopsIterations = estimations.size();
    double EstimationVal = 0;

    long double startEstimation =  (double) this->possibleGraphImagesOfPatternRoot.size();

    this->hopsEvaluation.hopsStd *= (startEstimation * startEstimation) / (long double) this->hopsEvaluation.hopsIterations;

    for (UInt64 Value : estimations){
        //Recalculate estimation from log value
        if (Value != std::numeric_limits<UInt64>::max()) {
            EstimationVal = exp((double) Value / (double) this->logValueMultiplier);
        }
            //Summarize all values and count zeros
        else{
            EstimationVal = 0;
            ++this->hopsEvaluation.hopsZeroIterations;
        }

        this->hopsEvaluation.hopsEstimation += EstimationVal;

        //Look at the different values of estimations
        if (!runParameters.single_number) {
            if (this->hopsEvaluation.estimationMap.find(std::round(EstimationVal) * startEstimation) ==
                this->hopsEvaluation.estimationMap.end()) {
                this->hopsEvaluation.estimationMap.insert({std::round(EstimationVal) * startEstimation, 1});
            } else {
                ++this->hopsEvaluation.estimationMap[std::round(EstimationVal) * startEstimation];
            }
        }
    }
    this->hopsEvaluation.hopsEstimation *= startEstimation / (long double) this->hopsEvaluation.hopsIterations;
    for (long double & snapshot : this->hopsEvaluation.snapshots) {
        snapshot *= startEstimation;
    }
    //Evaluation Runtime
    this->hopsEvaluation.evaluationRuntime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
    //Print and Save the results
    if (this->runParameters.print){
        this->hopsEvaluation.print();
    }
    if(runParameters.save && !this->out_path.empty()){
        this->hopsEvaluation.save(this->out_path);
    }
}
inline void Hops::EvaluateResult(int approximatedGED, HOPS_TYPE hopsType) {
    //Print and Save the results
    if (this->runParameters.print){
        this->hopsEvaluation.hopsEstimation = approximatedGED;
        this->hopsEvaluation.print(hopsType);
    }
    if(!this->out_path.empty()){
        this->hopsEvaluation.save(this->out_path);
    }
}

inline void Hops::InitEstimation(unsigned int& id, Hops::RunProps &runProps, const std::vector<NodeId>& PatternRootImages, bool random) {
    //Check if Id has reached integer limits, if so reset Id and all values for visited neighbors to 0
    if (id == std::numeric_limits<unsigned int>::max()){
        id = 0;
        for(NodeId& nodeId : runProps.savedNeighbors){
            nodeId = 0;
        }
    }
    ++id;
    if (random) {
        //Determine the initial graph image of the pattern root node
        NodeId rand = std::uniform_int_distribution<NodeId>(0, (NodeId) PatternRootImages.size() - 1)(runProps.gen);
        //Assign Initial Image
        runProps.currentRootNode = PatternRootImages[rand];
    }
    runProps.treeGraphMap[0] = runProps.currentRootNode;



    if (runProps.runType == RUN_TYPE::UNORDERED_SET) {
        runProps.graphImages.clear();
        runProps.graphImages.emplace(runProps.treeGraphMap[0]);
    }
    else if(runProps.runType == RUN_TYPE::VECTOR){
        runProps.savedNeighbors[runProps.treeGraphMap[0]] = id;
    }
}


inline bool Hops::UnlabeledEmbedding(unsigned int id, long double& accumulatedMean, long double& accumulatedDeviation, std::vector<UInt64>& estimationVector, UInt64& zeroIterations, RunProps& runProps) {
    NodeId CurrentGraphImageId;
    unsigned int SourceSize, TargetSize, possibleNeighbors;
    UInt64 estimation = 0;
    bool rootNode = true;
    // Iterate over all pairs of parent and child in the BFSOrder
    for (auto const & [parentNode, childrenNodes] : this->rootedPattern.BFSOrder) {
        CurrentGraphImageId = runProps.treeGraphMap[this->rootedPattern.GetBFSOrderIndexByNodeId(parentNode)];
        const Nodes& neighbors = this->currentGraph->get_neighbors(CurrentGraphImageId);
        SourceSize = static_cast<unsigned int>(childrenNodes.size());
        TargetSize = static_cast<unsigned int>(neighbors.size());
        if (rootNode) {
            if (SourceSize > TargetSize) {
                if (!this->runParameters.single_number) {
                    estimationVector.push_back(std::numeric_limits<UInt64>::max());
                }
                ++zeroIterations;
                return false;
            }
            rootNode = false;
        }
        else{
            if (SourceSize >= TargetSize) {
                if (!this->runParameters.single_number) {
                    estimationVector.push_back(std::numeric_limits<UInt64>::max());
                }
                ++zeroIterations;
                return false;
            }
        }

        //Check if there is valid neighbor mapping
        if (!Hops::PickRandomNeighbors(id, SourceSize, TargetSize, CurrentGraphImageId, neighbors, possibleNeighbors,
                                       childrenNodes, runProps)) {
            if (!this->runParameters.single_number) {
                estimationVector.push_back(std::numeric_limits<UInt64>::max());
            }
            ++zeroIterations;
            return false;
        }

        //Calculate the Value for estimation Combinations in Log-Space
        estimation += (this->preComputedLogValues[possibleNeighbors] - this->preComputedLogValues[possibleNeighbors - SourceSize]);

        if(!Hops::CheckValidGraphEmbedding(childrenNodes, runProps)){
            if (!this->runParameters.single_number) {
                estimationVector.push_back(std::numeric_limits<UInt64>::max());
            }
            ++zeroIterations;
            return false;
        }
    }
    if (this->runParameters.single_number) {
        long double val = exp((double) estimation / (double) this->logValueMultiplier);
        accumulatedMean += val;
        accumulatedDeviation += val*val;
    }
    else{
        estimationVector.push_back(estimation);
    }
    return true;
}

inline bool Hops::UnlabeledGraphEditDistance(unsigned int id, unsigned int &currentEditDistance, Hops::RunProps &runProps) {
    NodeId CurrentGraphImageId;
    unsigned int SourceSize, TargetSize;
    auto Error = static_cast<unsigned int>(this->currentPattern->nodes() + this->currentPattern->edges());
    unsigned int embeddedNodes = 1;
    unsigned int embeddedEdges = 0;
    for (auto const & [parentNode, childrenNodes] : this->rootedPattern.BFSOrder) {
        CurrentGraphImageId = runProps.treeGraphMap[this->rootedPattern.GetBFSOrderIndexByNodeId(parentNode)];
        SourceSize = static_cast<unsigned int>(childrenNodes.size());
        TargetSize = static_cast<unsigned int>(this->currentGraph->degree(CurrentGraphImageId));

        //Check if there is valid neighbor mapping
        if (!Hops::FindEmbedding(id, SourceSize, TargetSize, CurrentGraphImageId, embeddedNodes, embeddedEdges,
                                 childrenNodes, runProps)) {
            Hops::CheckValidGraphEmbedding(childrenNodes, embeddedEdges, runProps);
            Error -= (embeddedNodes+embeddedEdges);
            currentEditDistance = std::min(currentEditDistance, Error);
            return false;
        }

        if(!Hops::CheckValidGraphEmbedding(childrenNodes, embeddedEdges, runProps)){
            Error -= (embeddedNodes+embeddedEdges);
            currentEditDistance = std::min(currentEditDistance, Error);
            return false;
        }
    }
    currentEditDistance = 0;
    return true;
}


inline bool Hops::LabeledEmbeddings(unsigned int id, long double& accumulatedMean, long double& accumulatedDeviation, std::vector<UInt64>& estimationVector, UInt64& zeroIterations, RunProps& runProps) {
    NodeId PatternParentNode, GraphImageOfPatternParentNode;
    unsigned int PossibleNeighborsSize, SourceSize, TargetSize;
    UInt64 estimation = 0;
    bool rootNode = true;
    //Iterate over all nodes in the pattern graph
    for (int i = 0; i < this->rootedPattern.BFSOrder.size(); ++i) {
        PatternParentNode = this->rootedPattern.BFSOrder[i].first;
        GraphImageOfPatternParentNode = runProps.treeGraphMap[this->rootedPattern.GetBFSOrderIndexByNodeId(PatternParentNode)];
        const Nodes & neighbors = this->currentGraph->neighbors(GraphImageOfPatternParentNode);
        if(!this->currentPattern->isBigger(this->currentGraph->degreeByLabel(GraphImageOfPatternParentNode), PatternParentNode)) {
            for (auto const&[CurrentLabel, currentNodes] : this->rootedPattern.BFSOrderLabelMap[i]) {
                SourceSize = static_cast<int>(currentNodes.size());
                if (runProps.labelType == LABEL_TYPE::LABELED_DENSE) {
                    TargetSize = (!this->currentGraph->has_neighbor_label(GraphImageOfPatternParentNode, CurrentLabel)) ? 0
                                                                                                                        : static_cast<int>(this->currentGraph->degreeByLabel(GraphImageOfPatternParentNode)[CurrentLabel]);
                } else if (runProps.labelType == LABEL_TYPE::LABELED_SPARSE) {
                    TargetSize = static_cast<int>(this->currentGraph->degreeByLabel(GraphImageOfPatternParentNode)[CurrentLabel]);
                }
                if (SourceSize > TargetSize) {
                    if (!this->runParameters.single_number) {
                        estimationVector.push_back(std::numeric_limits<INDEX>::max());
                    }
                    ++zeroIterations;
                    return false;
                }

                //Do the embedding
                if(!PickRandomNeighbors(id, SourceSize, TargetSize,
                                        GraphImageOfPatternParentNode, neighbors,
                                        PossibleNeighborsSize, currentNodes, runProps,
                                        static_cast<int>(CurrentLabel))){
                    if (!this->runParameters.single_number) {
                        estimationVector.push_back(std::numeric_limits<INDEX>::max());
                    }
                    ++zeroIterations;
                    return false;
                }

                //Calculate the Value for estimation Combinations in Log-Space
                estimation += (this->preComputedLogValues[PossibleNeighborsSize] - this->preComputedLogValues[PossibleNeighborsSize - SourceSize]);
            }
        }
        else{
            if (!this->runParameters.single_number) {
                estimationVector.push_back(std::numeric_limits<INDEX>::max());
            }
            ++zeroIterations;
            return false;
        }

        //Check if all graph edges are valid
        if(!Hops::CheckValidGraphEmbedding(this->rootedPattern.BFSOrder[i].second, runProps)){
            if (!this->runParameters.single_number) {
                estimationVector.push_back(std::numeric_limits<INDEX>::max());
            }
            ++zeroIterations;
            return false;
        }
    }
    if (this->runParameters.single_number) {
        long double val = exp((double) estimation / (double) this->logValueMultiplier);
        accumulatedMean += val;
        accumulatedDeviation += val*val;
    }
    else{
        estimationVector.push_back(estimation);
    }
    return true;
}


inline bool Hops::LabeledGraphEditDistance(unsigned int id, unsigned int &currentEditDistance, Hops::RunProps &runProps) {
    NodeId PatternParentNode, GraphImageOfPatternParentNode;
    unsigned int SourceSize, TargetSize, PossibleNeighborsSize;
    bool rootNode = true;
    auto Error = static_cast<unsigned int>(this->currentPattern->nodes() + this->currentPattern->edges());
    unsigned int embeddedNodes = 1;
    unsigned int embeddedEdges = 0;
    //Iterate over all nodes in the pattern graph
    for (unsigned int i = 0; i < this->rootedPattern.BFSOrder.size(); ++i) {
        PatternParentNode = this->rootedPattern.BFSOrder[i].first;
        GraphImageOfPatternParentNode = runProps.treeGraphMap[this->rootedPattern.GetBFSOrderIndexByNodeId(PatternParentNode)];
        for (auto const&[CurrentLabel, currentNodes] : this->rootedPattern.BFSOrderLabelMap[i]) {
            SourceSize = static_cast<int>(currentNodes.size());
            if (runProps.labelType == LABEL_TYPE::LABELED_DENSE) {
                TargetSize = (!this->currentGraph->has_neighbor_label(GraphImageOfPatternParentNode, CurrentLabel))
                             ? 0
                             : static_cast<int>(this->currentGraph->degreeByLabel(GraphImageOfPatternParentNode)[CurrentLabel]);
            } else if (runProps.labelType == LABEL_TYPE::LABELED_SPARSE) {
                TargetSize = static_cast<int>(this->currentGraph->degreeByLabel(
                        GraphImageOfPatternParentNode)[CurrentLabel]);
            }
            //Do the embedding
            //Check if there is valid neighbor mapping
            if (!Hops::FindEmbedding(id, SourceSize, TargetSize, GraphImageOfPatternParentNode, embeddedNodes,
                                     embeddedEdges,
                                     currentNodes, runProps, static_cast<int>(CurrentLabel))) {
                Hops::CheckValidGraphEmbedding(currentNodes, embeddedEdges, runProps);
                Error -= (embeddedNodes + embeddedEdges);
                currentEditDistance = std::min(currentEditDistance, Error);
                return false;
            }
        }

        //Check if all graph edges are valid
        if(!Hops::CheckValidGraphEmbedding(this->rootedPattern.BFSOrder[i].second, runProps)){
            Error -= (embeddedNodes + embeddedEdges);
            currentEditDistance = std::min(currentEditDistance, Error);
            return false;
        }
    }
    currentEditDistance = 0;
    return true;
}


inline bool Hops::PickRandomNeighbors(unsigned int id, unsigned int sourceSize, unsigned int targetSize, NodeId currentNodeId, const Nodes& currentNodeNeighbors, unsigned int& possibleNeighbors, const  Nodes &orderNodes, RunProps& runProps, int label) const {
    possibleNeighbors = 0;
    unsigned int rand_number, Successes = 0;
    NodeId graphImageNode;
    long int allowedFailures = targetSize - sourceSize;
    unsigned int maxBound = targetSize - 1;
    for (int i = 0; i < targetSize; ++i) {
        if (Successes < sourceSize) {
            rand_number = std::uniform_int_distribution<unsigned int>(Successes, maxBound)(runProps.gen);
            if (runProps.labelType != LABEL_TYPE::UNLABELED) {
                graphImageNode = this->currentGraph->neighbor(currentNodeId, runProps.randomNeighbors[rand_number], label);
            }
            else{
                graphImageNode = currentNodeNeighbors[runProps.randomNeighbors[rand_number]];
            }
            if (runProps.VisitNeighbor(id, graphImageNode)) {
                //Random neighbor is not mapped yet
                runProps.swapNeighbors(rand_number, Successes);
                NodeId patternVertex = this->rootedPattern.GetBFSOrderIndexByNodeId(orderNodes[Successes]);
                if (this->currentPattern->degree(orderNodes[Successes]) > this->currentGraph->degree(graphImageNode)){
                    //Swap back vector for random neighbor assignment
                    runProps.swapBackNeighbors();
                    return false;
                }
                else{
                    //Map patternVertex to graphImageNode
                    runProps.treeGraphMap[patternVertex] = graphImageNode;
                    ++Successes, ++possibleNeighbors;
                }
            } else {
                //Random neighbor is already mapped
                --allowedFailures;
                if (allowedFailures < 0) {
                    //Swap back vector for random neighbor assignment
                    runProps.swapBackNeighbors();
                    return false;
                }
                runProps.swapNeighbors(rand_number, maxBound);
                --maxBound;
            }
        }
        else {
            break;
        }
    }
    //Check if remaining neighbors could have been possible neighbors
    for (unsigned int i = Successes; i <= maxBound; ++i) {
        if (runProps.labelType != LABEL_TYPE::UNLABELED) {
            graphImageNode = this->currentGraph->neighbor(currentNodeId, runProps.randomNeighbors[i],
                                                          label);
        }
        else{
            graphImageNode = this->currentGraph->neighbor(currentNodeId, runProps.randomNeighbors[i]);
        }
        if(!runProps.IsNeighborVisited(id, graphImageNode)){
            ++possibleNeighbors;
        }
    }
    //Swap back vector for random neighbor assignment
    runProps.swapBackNeighbors();
    return Successes == sourceSize;
}

inline bool Hops::FindEmbedding(unsigned int id, unsigned int sourceSize, unsigned int targetSize, NodeId currentNodeId, unsigned int &numberOfEmbeddedNodes, unsigned int &numberOfEmbeddedEdges, const Nodes &orderNodes, Hops::RunProps &runProps,
                         int label) const {
    unsigned int Successes = 0;
    unsigned int rand_number;
    NodeId NeighborNodeId;
    unsigned int maxBound = targetSize - 1;
    bool success = true;
    for (int i = 0; i < targetSize; ++i) {
        if (Successes < sourceSize) {
            rand_number = std::uniform_int_distribution<unsigned int>(Successes, maxBound)(runProps.gen);
            if (runProps.labelType != LABEL_TYPE::UNLABELED) {
                NeighborNodeId = this->currentGraph->neighbor(currentNodeId, runProps.randomNeighbors[rand_number],
                                                              label);
            }
            else{
                NeighborNodeId = this->currentGraph->neighbor(currentNodeId, runProps.randomNeighbors[rand_number]);
            }
            if (runProps.VisitNeighbor(id, NeighborNodeId)) {
                //Random neighbor is not mapped yet
                runProps.swapNeighbors(rand_number, Successes);
                runProps.treeGraphMap[this->rootedPattern.GetBFSOrderIndexByNodeId(orderNodes[Successes])] = NeighborNodeId;
                ++Successes, ++numberOfEmbeddedNodes, ++numberOfEmbeddedEdges;
            } else {
                //Random neighbor is already mapped
                runProps.swapNeighbors(rand_number, maxBound);
                --maxBound;
            }
        }
        else{
            break;
        }
    }
    //Swap back vector for random neighbor assignment
    runProps.swapBackNeighbors();
    return Successes == sourceSize;
}


inline bool Hops::CheckValidGraphEmbedding(const Nodes &nodes, const RunProps& runProps) const {
    NodeId sourceId, targetId;
    //Check if all graph edges are valid
    if (!this->currentPatternIsTree) {
        for (NodeId nodeId: nodes) {
            for (auto &[source, target]: this->nonTreeEdges[this->rootedPattern.GetBFSOrderIndexByNodeId(nodeId)]) {
                sourceId = runProps.treeGraphMap[source];
                targetId = runProps.treeGraphMap[target];

                const Nodes &sourceNeighbors = this->currentGraph->neighbors(sourceId);
                const Nodes &targetNeighbors = this->currentGraph->neighbors(targetId);

                if (sourceNeighbors.size() <= targetNeighbors.size()) {
//                    if(!std::binary_search(sourceNeighbors.begin(), sourceNeighbors.end(), targetId)){
//                        return false;
//                    }
                    if (std::find(sourceNeighbors.begin(), sourceNeighbors.end(), targetId) == sourceNeighbors.end()) {
                        return false;
                    }
                } else {
//                    if(!std::binary_search(targetNeighbors.begin(), targetNeighbors.end(), sourceId)){
//                        return false;
//                    }
                    if (std::find(targetNeighbors.begin(), targetNeighbors.end(), sourceId) == targetNeighbors.end()) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

inline bool Hops::CheckValidGraphEmbedding(const Nodes &nodes, unsigned int &embeddedEdges, const Hops::RunProps &runProps) const {
    NodeId sourceId, targetId;
    //Check if all graph edges are valid
    if (!this->currentPatternIsTree)
    {
        int EdgeError = 0;
        for (NodeId nodeId : nodes)
        {
            for (auto & [source, target] : this->nonTreeEdges[this->rootedPattern.GetBFSOrderIndexByNodeId(nodeId)]) {
                sourceId = runProps.treeGraphMap[source];
                targetId = runProps.treeGraphMap[target];
                const Nodes& sourceNeighbors = this->currentGraph->neighbors(sourceId);
                const Nodes& targetNeighbors = this->currentGraph->neighbors(targetId);
                if (sourceNeighbors.size() < targetNeighbors.size()){
                    if(std::binary_search(sourceNeighbors.begin(), sourceNeighbors.end(), targetId)){
                        ++embeddedEdges;
                    }
                    else{
                        ++EdgeError;
                    }
                }
                else{
                    if(std::binary_search(targetNeighbors.begin(), targetNeighbors.end(), sourceId)){
                        ++embeddedEdges;
                    }
                    else{
                        ++EdgeError;
                    }
                }
            }
        }
        return EdgeError == 0;
    }
    return true;
}

void Hops::Automorphisms(size_t graphId,RunParameters rParameters) {
    GraphStruct pattern = this->graphs[graphId];
    this->Run(graphId, pattern, rParameters);
}


inline bool Hops::RunProps::IsNeighborVisited(unsigned int Id, NodeId neighborNodeId) {
    if (runType == RUN_TYPE::VECTOR){
        return this->savedNeighbors[neighborNodeId] == Id;
    }
    else{
        return this->graphImages.find(neighborNodeId) != this->graphImages.end();
    }
}

inline bool Hops::RunProps::VisitNeighbor(unsigned int Id, NodeId neighborNodeId) {
    if (runType == RUN_TYPE::VECTOR){
        if(this->savedNeighbors[neighborNodeId] != Id){
            this->savedNeighbors[neighborNodeId] = Id;
            return true;
        }
        return false;
    }
    else{
        if(this->graphImages.find(neighborNodeId) == this->graphImages.end()){
            this->graphImages.emplace(neighborNodeId);
            return true;
        }
        return false;
    }
}

inline void Hops::RunProps::swapNeighbors(unsigned int idx_a, unsigned int idx_b) {
    std::swap(randomNeighbors[idx_a], randomNeighbors[idx_b]);
    randomNeighborsSwapPairs[swapCount].first=idx_a;
    randomNeighborsSwapPairs[swapCount].second=idx_b;
    ++swapCount;
}

inline void Hops::RunProps::swapBackNeighbors() {
    for (unsigned int i = swapCount; i > 0; --i) {
        const auto& [n,m] = randomNeighborsSwapPairs[i-1];
        std::swap(randomNeighbors[n], randomNeighbors[m]);
    }
    swapCount = 0;
}


#endif //HOPS_HOPS_H
