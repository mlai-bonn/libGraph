//
// Created by Florian on 12.04.2021.
//

#ifndef CLOSURES_GRAPHCLOSURES_H
#define CLOSURES_GRAPHCLOSURES_H


#include <string>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <unordered_set>
#include "BaseOperator.h"

enum class EGraphClosureType {
    EXACT_GEODESIC,
    APPROXIMATION_PRECLOSURE,
    APPROXIMATION_OUTERPLANAR,
    APPROXIMATION_TREE,
    APPROXIMATION_LAYERING,
};

/**
 * @brief Parameters for the closure computation on graphs
 * @details The parameters are used to define the closure computation
 * @see ClosureParameters
 * @param closureType defines the closure function
 * @param number_of_substructures defines the number of substructures used for the approximation
 * @param relative_threshold defines the relative threshold for the approximation
 */
struct GraphClosureParameters : ClosureParameters {
    EGraphClosureType closureType = EGraphClosureType::EXACT_GEODESIC;

    // Used only for approximation methods
    int number_of_substructures = 1;
    std::string path_to_substructures;
    double relative_threshold = 0.0;
};

/**
 * @brief Class for computing the closure of a _graph
 */
class GraphClosure : public BaseOperator<GraphStruct> {
public:
    /**
     * @brief Default constructor
     * @param graph Input _graph
     */
    explicit GraphClosure(GraphStruct& graph) : _graph(graph){
        Id = 0;
        _graph_distances = std::vector<int>(graph.nodes(), -2);
        _graph_containment_list = std::vector<int>(graph.nodes(), 0);
        _graph_predecessors = std::vector<std::vector<NodeId>>(graph.nodes(), std::vector<NodeId>());
        _graph_bfs_list = std::vector<int>(graph.nodes(), 0);
    };

    /**
     * @brief Computes the closure of the input set, the closure function is determined in the closure parameters
     * @param closureParameters defines all necessary parameters for the closure computation
     */
    void closure(ClosureParameters& closureParameters) override;
    void exact_geodesic_closure(GraphClosureParameters& closureParameters);

    /**
     * @brief Computes the closure of the input set
     * @param outerplanarGraphData
     * @param closureParameters
     */
    void outerplanar_closure(GraphClosureParameters &closureParameters);


private:
    /**
     * @brief Forward search through the _graph
     * @param graph the _graph
     * @param closureParameters defines all necessary parameters for the closure computation
     */
    void bfs_forward(GraphClosureParameters& closureParameters, NodeId bfs_start, std::deque<NodeId> &bfsQueue);
    /**
     * @brief Backward search
     * @param closureParameters
     * @param bfsQueue
     * @param bfsElements
     * @param generatedElements
     */
    void bfs_backward(GraphClosureParameters &closureParameters,
                      std::deque<NodeId> &bfsQueue, std::deque<NodeId>& bfsElements, std::unordered_set<NodeId>& generatedElements);

    void set_break_condition(bool &break_condition, std::deque<NodeId> &newElements, GraphClosureParameters &closureParameters);

    static void get_generators(OuterplanarComponent& outerplanarComponent, std::set<NodeId> &input_set,
                        std::set<NodeId> &generator_set);

private:
    int Id = 0;
    std::vector<int> _graph_distances;
    std::vector<std::vector<NodeId>> _graph_predecessors;
    std::vector<int> _graph_containment_list;
    std::vector<int> _graph_bfs_list;
    GraphStruct& _graph;
    std::chrono::time_point<std::chrono::system_clock> time;
    int _iterations = 0;
};

inline void GraphClosure::closure(ClosureParameters& closureParameters) {
    // check if closureParameters can be casted to GraphClosureParameters
    auto* graphClosureParameters = (GraphClosureParameters*)(&closureParameters);
    if (graphClosureParameters != nullptr){
        switch (graphClosureParameters->closureType) {
            case EGraphClosureType::EXACT_GEODESIC:
                exact_geodesic_closure(*graphClosureParameters);
                break;
            case EGraphClosureType::APPROXIMATION_PRECLOSURE:
                graphClosureParameters->onlyPreClosure = true;
                exact_geodesic_closure(*graphClosureParameters);
                break;
            case EGraphClosureType::APPROXIMATION_OUTERPLANAR:
                break;
            case EGraphClosureType::APPROXIMATION_TREE:
                break;
            case EGraphClosureType::APPROXIMATION_LAYERING:
                break;
        }
    }
    else{
        // throw an error message that one should use GraphClosureParameters instead of ClosureParameters


    }
}

inline void GraphClosure::exact_geodesic_closure(GraphClosureParameters& closureParameters) {
    // Check validity of input set
    if (closureParameters.input_set.empty()) {
        return;
    }

    //Check if outerplanar closure should be computed
    if (dynamic_cast<OuterplanarGraphData *>(&_graph) != nullptr) {
        outerplanar_closure(closureParameters);
        return;
    }

    // If the _graph type is tree or outerplanar only preclosure is needed
    if (_graph.GetType() == GraphType::TREE || _graph.GetType() == GraphType::OUTERPLANAR) {
        closureParameters.onlyPreClosure = true;
    }

    //Initialize variables
    _graph_distances.resize(_graph.nodes(), -2);
    std::fill(_graph_distances.begin(), _graph_distances.end(), -2);
    _graph_containment_list.resize(_graph.nodes(), 0);
    _graph_predecessors.resize(_graph.nodes(), std::vector<NodeId>());
    std::fill(_graph_predecessors.begin(), _graph_predecessors.end(), std::vector<NodeId>());
    _graph_bfs_list.resize(_graph.nodes(), 0);


    std::unordered_set<NodeId> generatedElements;
    std::deque<NodeId> bfsQueue;
    // set closed set to input set
    closureParameters.closed_set.clear();
    closureParameters.closed_set = closureParameters.input_set;
    // clear added elements
    closureParameters.added_elements.clear();
    //Clear params
    //closureParameters.preClosureDepthToElements.clear();

    //Variable initialization
    std::deque<NodeId> newElements = std::deque<NodeId>(closureParameters.input_set.begin(),
                                                        closureParameters.input_set.end());

    closureParameters.preclosure_depth = 0;
    NodeId LastElemOfPreviousStep = newElements.back();
    _iterations = 0;
    bool break_condition = false;
    this->time = std::chrono::system_clock::now();
    set_break_condition(break_condition, newElements, closureParameters);
    closureParameters.preclosure_depth = 1;
    NodeId last_pre_closure_element = newElements.back();
    bool new_pre_closure_step = false;
    while (!break_condition) {
        NodeId bfsStart = newElements.front();

//Set detailed closure analysis
//            if (closureParameters.detailed_analysis) {
//                if (closureParameters.preClosureDepthToElements.find(closureParameters.preclosure_depth) ==
//                    closureParameters.preClosureDepthToElements.end()) {
//                    // insert preclosure depth and empty vector to the map preClosureDepthToElements
//                    auto insert_pair = std::pair<int, std::vector<int>>(closureParameters.preclosure_depth,
//                                                                        std::vector<int>());
//                    //closureParameters.preClosureDepthToElements.insert(insert_pair);
//                }
//                closureParameters.preClosureDepthToElements[closureParameters.preclosure_depth].emplace_back(bfsStart);
//            }
        if (bfsStart == last_pre_closure_element) {
            new_pre_closure_step = true;
        }
        if (new_pre_closure_step) {
            last_pre_closure_element = newElements.back();
            new_pre_closure_step = false;
            ++closureParameters.preclosure_depth;
        }

        newElements.pop_front();
//                if (_graph.graphType == graphType::OUTERPLANAR) {
//                    auto const it = generatedElements.find(bfsStart);
//                    if (it == generatedElements.end()) {
//                        bfs_forward(_graph, closureParameters.closed_set, bfsStart, bfsQueue);
//                        bfs_backward(_graph, closureParameters.closed_set, closureParameters, bfsQueue, newElements,
//                                     generatedElements);
//                        ++iterations;
//                    }
//                } else {
        bfs_forward(closureParameters, bfsStart, bfsQueue);
        bfs_backward(closureParameters, bfsQueue, newElements, generatedElements);
        ++_iterations;
//}
        // if the _graph is a tree and the threshold is infinite one bfs is enough to find the closure
        if (_graph.GetType() == GraphType::TREE && closureParameters.threshold == std::numeric_limits<int>::max()) {
            break;
        }
        set_break_condition(break_condition, newElements, closureParameters);
    }
}

inline void GraphClosure::bfs_forward(GraphClosureParameters& closureParameters, NodeId bfs_start, std::deque<NodeId> &bfsQueue){

    //Forward Search through _graph O(m)
    //all elements in closed set
    ++Id;
    if (Id == 0){
        std::fill(_graph_containment_list.begin(), _graph_containment_list.end(), 0);
        std::fill(_graph_bfs_list.begin(), _graph_bfs_list.end(), 0);
        ++Id;
    }
    for (NodeId elem : closureParameters.closed_set) {
        _graph_containment_list[elem] = Id;
    }
    _graph_distances[bfs_start] = 0;
    _graph_bfs_list[bfs_start] = Id;
    _graph_predecessors[bfs_start].clear();
    bfsQueue.clear();
    bfsQueue.push_back(bfs_start);
    int visitedSize = 1;
    int lastDistance = std::numeric_limits<int>::max();
    while (!bfsQueue.empty()) {
        NodeId currentNodeId = bfsQueue.back();
        bfsQueue.pop_back();
        int currentDistance = _graph_distances[currentNodeId];
        if ((visitedSize == closureParameters.closed_set.size() && _graph_distances[currentNodeId] > lastDistance)) {
            break;
        }
        // Use threshold to stop search this corresponds to considering only shortest paths of length <= threshold
        if (currentDistance + 1 > closureParameters.threshold){
            break;
        }
        for (auto neighborId : _graph.get_neighbors(currentNodeId)) {
            //Neighbor is unvisited by bfs
            if (_graph_bfs_list[neighborId] != Id) {
                _graph_distances[neighborId] = currentDistance + 1;
                _graph_bfs_list[neighborId] = Id;
                _graph_predecessors[neighborId].clear();
                bfsQueue.push_front(neighborId);
                //Neighbor is in closed set
                if (_graph_containment_list[neighborId] == Id) {
                    ++visitedSize;
                }
            }
            if (_graph_distances[neighborId] == currentDistance + 1) {
                _graph_predecessors[neighborId].push_back(currentNodeId);
            }
        }
        if (visitedSize == closureParameters.closed_set.size()) {
            lastDistance = currentDistance;
        }
    }
}

inline void GraphClosure::bfs_backward(GraphClosureParameters &closureParameters,
                                       std::deque<NodeId> &bfsQueue, std::deque<NodeId>& bfsElements, std::unordered_set<NodeId>& generatedElements) {
    //Backward search
    ++Id;
    if (Id == 0){
        std::fill(_graph_containment_list.begin(), _graph_containment_list.end(), 0);
        ++Id;
    }
    bfsQueue.clear();
//        if (_graph.graphType == GraphType::OUTERPLANAR){
//            for (NodeId elem: closureParameters.input_set) {
//                _graph_containment_list[elem] = _Id;
//                bfsQueue.push_back(elem);
//            }
//        }
//        else {
    for (NodeId elem: closureParameters.closed_set) {
        _graph_containment_list[elem] = Id;
        bfsQueue.push_back(elem);
    }
//        }

    while (!bfsQueue.empty()){
        NodeId NextElementId = bfsQueue.back();
        bfsQueue.pop_back();
        for (NodeId NeighborNodeId : _graph_predecessors[NextElementId]) {
            if (_graph_containment_list[NeighborNodeId] != Id) {
                _graph_containment_list[NeighborNodeId] = Id;
                closureParameters.closed_set.insert(NeighborNodeId);
                closureParameters.added_elements.insert(NeighborNodeId);
                bfsQueue.push_back(NeighborNodeId);

                // Add the newly found elements to the bfs start nodes if not only the preclosure is computed
                if (!closureParameters.onlyPreClosure){
                    bfsElements.push_back(NeighborNodeId);
                }
            }
            else{
                if (_graph.GetType() == GraphType::TREE){
                    generatedElements.insert(NeighborNodeId);
                }
            }
        }
    }

}

inline void GraphClosure::outerplanar_closure(GraphClosureParameters &closureParameters) {
    auto* graph = dynamic_cast<OuterplanarGraphData *>(&_graph);
    if (graph != nullptr) {

        GraphClosureParameters bbTreeParameters;
        for (const auto &elem: closureParameters.input_set) {
            graph->get_bb_tree_ids(elem, bbTreeParameters.input_set);
        }
        GraphClosure closure_bb_tree = GraphClosure(graph->get_bbTree());
        closure_bb_tree.closure(bbTreeParameters);
        int nodeId = 0;
        int componentId = 0;
        std::vector<std::set<NodeId>> componentInput = std::vector<std::set<NodeId>>(graph->Components.size(),
                                                                               std::set<NodeId>());
        std::vector<NodeOrComponent> nodeOrComponents;
        for (auto const &inputElement: closureParameters.input_set) {
            graph->get_components(inputElement, nodeOrComponents);
            for (auto const &nodeOrComp: nodeOrComponents) {
                if (nodeOrComp.is_component(componentId)) {
                    OuterplanarComponent &currentComponent = graph->Components[componentId];
                    componentInput[componentId].insert( currentComponent.NodeIdToComponentNodeId[inputElement]);
                }
            }
        }
        for (auto const &elem: bbTreeParameters.closed_set) {
            NodeOrComponent nodeComponent = graph->get_bbNodeOrComponent(elem);
            if (nodeComponent.is_node(nodeId)) {
                closureParameters.closed_set.insert(nodeId);
                graph->get_components(nodeId, nodeOrComponents);
                for (auto const &nodeOrComp: nodeOrComponents) {
                    if (nodeOrComp.is_component(componentId)) {
                        OuterplanarComponent &currentComponent = graph->Components[componentId];
                        componentInput[componentId].insert((int) currentComponent.componentId(nodeId));
                    }
                }
            }
        }

        for (int i = 0; i < graph->Components.size(); ++i) {
            OuterplanarComponent &component = graph->Components[i];
            GraphClosureParameters componentGraphClosureParameters;
            if (componentInput[i].size() > 1) {
                GraphClosure::get_generators(component, componentInput[i], componentGraphClosureParameters.input_set);
                GraphClosure componentClosure = GraphClosure(component.component);
                componentClosure.closure(componentGraphClosureParameters);
                for (auto const &elem: componentGraphClosureParameters.closed_set) {
                    closureParameters.closed_set.insert(component.nodeId(elem));
                }
            }
        }
    }
}

inline void GraphClosure::set_break_condition(bool &break_condition, std::deque<NodeId> &newElements, GraphClosureParameters &closureParameters)
{
    break_condition = newElements.empty();
    if (closureParameters.iteration_number > 0){
        break_condition = break_condition || _iterations > closureParameters.iteration_number;
    }
    if (closureParameters.timeConstraint != -1){
        double time_count = (double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - this->time).count();
        double max_time =  closureParameters.timeConstraint * 1000000;
        break_condition = break_condition || time_count > max_time ;
    }
}

void GraphClosure::get_generators(OuterplanarComponent& outerplanarComponent, std::set<NodeId> &input_set,
                                  std::set<NodeId> &generator_set) {
    // TODO get the three defining generators for each face
    generator_set.clear();
    std::vector<GraphStruct> faces;
    std::vector<std::vector<NodeId>> nodeToFaces;
    generator_set = input_set;
}



#endif //CLOSURES_GRAPHCLOSURES_H
