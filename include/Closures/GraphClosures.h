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
#include "../../include/io/StaticFunctions.h"

typedef void (*ClosureFunction)(GraphStruct&, ClosureParameters&);

class GraphClosureSP : public BaseOperator<GraphStruct> {
public:
    explicit GraphClosureSP(const GraphStruct& graph) : _graph(graph){
        _Id = 0;
        _graph_distances = std::vector<int>(graph.nodes(), -2);
        _graph_containment_list = std::vector<int>(graph.nodes(), 0);
        _graph_predecessors = std::vector<std::vector<NodeId>>(graph.nodes(), std::vector<NodeId>());
        _graph_bfs_list = std::vector<int>(graph.nodes(), 0);
        _is_tree = graph.CheckTree();
    };

    void closure(ClosureParameters& closureParameters) override {
        if (closureParameters.input_set.empty()) {
            return;
        }
        _graph_distances.resize(_graph.nodes(), -2);
        _graph_containment_list.resize(_graph.nodes(), 0);
        _graph_predecessors.resize(_graph.nodes(), std::vector<NodeId>());
        _graph_bfs_list.resize(_graph.nodes(), 0);


        std::unordered_set<NodeId> generatedElements;
        std::deque<NodeId> bfsQueue;
        // set closed set to input set
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
        NodeId last_preclosure_element = newElements.back();
        bool new_preclosure_step = false;
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
            if (bfsStart == last_preclosure_element) {
                new_preclosure_step = true;
            }
            if (new_preclosure_step) {
                last_preclosure_element = newElements.back();
                new_preclosure_step = false;
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
            bfs_forward(closureParameters.closed_set, bfsStart, bfsQueue);
            bfs_backward(closureParameters.closed_set, closureParameters, bfsQueue, newElements, generatedElements);
            ++_iterations;
            //}
            if (_is_tree) {
                break;
            }
            set_break_condition(break_condition, newElements, closureParameters);
        }
    }


    std::pair<int, int> _info;
private:
    void bfs_forward(std::set<NodeId>& target_set, NodeId bfs_start,
                std::deque<NodeId> &bfsQueue) {

        //Forward Search through graph O(m)
        //all elements in closed set
        ++_Id;
        if (_Id == 0){
            std::fill(_graph_containment_list.begin(), _graph_containment_list.end(), 0);
            std::fill(_graph_bfs_list.begin(), _graph_bfs_list.end(), 0);
            ++_Id;
        }
        for (NodeId elem : target_set) {
            _graph_containment_list[elem] = _Id;
        }
        _graph_distances[bfs_start] = 0;
        _graph_bfs_list[bfs_start] = _Id;
        _graph_predecessors[bfs_start].clear();
        bfsQueue.clear();
        bfsQueue.push_back(bfs_start);
        int visitedSize = 1;
        int lastDistance = std::numeric_limits<int>::max();
        while (!bfsQueue.empty()) {
            NodeId currentNodeId = bfsQueue.back();
            if (_info.first == 17 && _info.second == 5 && _iterations == 26) {
                // print current node id
                std::cout << "Current Node Id: " << currentNodeId << std::endl;
            }
            bfsQueue.pop_back();
            int currentDistance = _graph_distances[currentNodeId];
            if ((visitedSize == target_set.size() && _graph_distances[currentNodeId] > lastDistance)) {
                break;
            }
            for (auto neighborId : _graph.get_neighbors(currentNodeId)) {
                if (_info.first == 17 && _info.second == 5 && _iterations == 26) {
                    std::cout << "Neighbor Id: " << currentNodeId << std::endl;
                }
                //Neighbor is unvisited by bfs
                if (_graph_bfs_list[neighborId] != _Id) {
                    _graph_distances[neighborId] = currentDistance + 1;
                    _graph_bfs_list[neighborId] = _Id;
                    _graph_predecessors[neighborId].clear();
                    bfsQueue.push_front(neighborId);
                    //Neighbor is in closed set
                    if (_graph_containment_list[neighborId] == _Id) {
                        ++visitedSize;
                    }
                }
                if (_graph_distances[neighborId] == currentDistance + 1) {
                    _graph_predecessors[neighborId].push_back(currentNodeId);
                }
            }
            if (visitedSize == target_set.size()) {
                lastDistance = currentDistance;
            }
        }
    }

    void bfs_backward(std::set<NodeId>& target_set, ClosureParameters &closureParameters,
                                      std::deque<NodeId> &bfsQueue, std::deque<NodeId>& bfsElements, std::unordered_set<NodeId>& generatedElements) {
        //Backward search
        ++_Id;
        if (_Id == 0){
            std::fill(_graph_containment_list.begin(), _graph_containment_list.end(), 0);
            ++_Id;
        }
        bfsQueue.clear();
//        if (graph.graphType == GraphType::OUTERPLANAR){
//            for (NodeId elem: closureParameters.input_set) {
//                _graph_containment_list[elem] = _Id;
//                bfsQueue.push_back(elem);
//            }
//        }
//        else {
        for (NodeId elem: target_set) {
            _graph_containment_list[elem] = _Id;
            bfsQueue.push_back(elem);
        }
//        }

        while (!bfsQueue.empty()){
            NodeId NextElementId = bfsQueue.back();
            bfsQueue.pop_back();
            for (NodeId NeighborNodeId : _graph_predecessors[NextElementId]) {
                if (_graph_containment_list[NeighborNodeId] != _Id) {
                    _graph_containment_list[NeighborNodeId] = _Id;
                    closureParameters.closed_set.insert(NeighborNodeId);
                    closureParameters.added_elements.insert(NeighborNodeId);
                    bfsQueue.push_back(NeighborNodeId);

                    //New elements should be only considered in the case that the graph is of general type
                    if (!_is_tree) {
                        bfsElements.push_back(NeighborNodeId);
                    }
                }
                else{
                    if (_is_tree){
                        generatedElements.insert(NeighborNodeId);
                    }
                }
            }
        }

    }


    void set_break_condition(bool &break_condition, std::deque<NodeId> &newElements, ClosureParameters &closureParameters)
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

private:
    int _Id = 0;
    std::vector<int> _graph_distances;
    std::vector<std::vector<NodeId>> _graph_predecessors;
    std::vector<int> _graph_containment_list;
    std::vector<int> _graph_bfs_list;
    const GraphStruct& _graph;
    int _threshold = std::numeric_limits<int>::max();
    std::chrono::time_point<std::chrono::system_clock> time;

    bool _is_tree = false;

    int _iterations = 0;
};


#endif //CLOSURES_GRAPHCLOSURES_H
