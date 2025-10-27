//
// Created by Florian on 09.03.2021.
//

#ifndef HOPS_GRAPHFUNCTIONS_H
#define HOPS_GRAPHFUNCTIONS_H

#include <random>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <queue>
#include "GraphDataStructures/GraphStructs.h"


class GraphFunctions {
public:
    static std::unordered_map<Label, INDEX> GetLabelFrequency(std::unordered_map<Label, Nodes>& labelMap);
    static std::unordered_map<Label, Nodes> GetGraphLabelMap(Labels& nodeLabels);
    static std::pair<std::unordered_map<Label, Label>, std::unordered_map<Label, Label>> GetLabelMappingPair(Labels& nodeLabels);
    static void GetNodesByLabel(Labels& nodeLabels, Label label, Nodes& nodesByLabel);
    static Labels GetLabelsByNodes(Nodes& nodeIds, const Labels& graphLabels);
    static std::unordered_set<Label> GetUniqueLabelsByNodes(Nodes& nodeIds, const Labels& GraphLabels);

    template<typename T>
    static void GetComponents(T &graph, std::vector<INDEX> &components);

    template<typename T>
    static void GetComponents(T &graph, std::vector<Nodes> &components);

    template<typename T>
    static void GetLargestComponent(T &graph, Nodes &nodes);

    template<typename T>
    static T SubGraph(const T &graph, const Nodes &nodeIds);

    template<typename T>
    static T GetLargestComponent(T &graph);

    template <typename T>
    static bool ReachableNodes(T& graph, NodeId root, std::vector<INDEX>& reachability, INDEX Id, INDEX& number);

    template<typename T>
    static void BFSDistances(const T &graph, INDEX root, std::vector<INDEX> &distances);
    template<typename T>
    static void DFS(const T &graph, T &tree, Nodes& nodeOrder, NodeId rootNodeId = -1, int seed = 0);
    template<typename T>
    static void OptOrdering(const T &graph, Nodes& nodeOrder);
    template<typename T>
    static void ReorderGraph(T &graph, const Nodes& nodeOrder);

    template<typename T1, typename T2>
    static bool CheckSpanningTree(const T1 &graph, T2 &spanningTree);
    template<typename T>
    static bool CheckTree(const T& graph);

    
};

template<typename T>
bool GraphFunctions::ReachableNodes(T &graph, const NodeId root, std::vector<INDEX>& reachability, const INDEX Id, INDEX& number){
    number = 1;
    reachability[root] = Id;
    std::deque<NodeId> nodes;
    nodes.push_back(root);
    while (!nodes.empty()){
        NodeId currentNode = nodes.back();
        nodes.pop_back();
        for (NodeId i = 0; i < graph.degree(currentNode); ++i) {
            NodeId neighbor = graph.neighbor(currentNode, i);
            if (reachability[neighbor] != Id){
                reachability[neighbor] = Id;
                nodes.push_front(neighbor);
                ++number;
            }
        }
    }
    return number == graph.nodes();
}

template<typename T>
 void GraphFunctions::GetComponents(T &graph, std::vector<INDEX> & components){
    components.clear();
    components.resize(graph.nodes(), -1);
    INDEX number = 0;
    INDEX pos = 0;
    INDEX Id = 0;
    while (number < graph.nodes()){
        INDEX compSize = 0;
        while(components[pos] != -1){
            ++pos;
        }
        ReachableNodes(graph, pos, components, Id, compSize);
        number += compSize;
        ++Id;
        ++pos;
    }

}

template<typename T>
 void GraphFunctions::GetComponents(T &graph, std::vector<Nodes> & components){
    std::vector<INDEX> com;
    GraphFunctions::GetComponents(graph, com);
    for (INDEX i = 0; i < com.size(); ++i) {
        INDEX componentId = com[i];
        if (components.empty() || components.size() <= componentId){
            components.resize(componentId + 1, Nodes());
        }
        components[componentId].emplace_back(i);
    }
}

template<typename T>
 void GraphFunctions::GetLargestComponent(T &graph, Nodes &nodes) {
    std::vector<Nodes> components;
    GetComponents(graph, components);
    INDEX maxSize = 0;
    INDEX maxIndex = 0;
    for (INDEX i = 0; i < components.size(); ++i) {
        if (const auto compSize = (INDEX) components[i].size(); compSize > maxSize){
            maxSize = compSize;
            maxIndex = i;
        }
    }
    nodes = components[maxIndex];
}

template<typename T>
 T GraphFunctions::SubGraph(const T &graph, const Nodes& nodeIds){
    std::string graph_name = graph.GetName() + "_subgraph";
    T g = T(graph_name, static_cast<INDEX>(nodeIds.size()));
    for (INDEX i = 0; i < nodeIds.size(); ++i) {
        g.IdsToOriginalIds.emplace(nodeIds[i], i);
    }
    for (INDEX srcNode : nodeIds) {
        for (INDEX dstNode : graph.graph()[srcNode]) {
            if (std::ranges::find(nodeIds, dstNode) != nodeIds.end()){
                g.AddEdge(g.IdsToOriginalIds[srcNode], g.IdsToOriginalIds[dstNode], true);
            }
        }
    }
    return g;
}

template<typename T>
 T GraphFunctions::GetLargestComponent(T &graph){
    Nodes nodes;
    GraphFunctions::GetLargestComponent(graph, nodes);
    return SubGraph(graph, nodes);
}

template<typename T>
inline void GraphFunctions::BFSDistances(const T &graph, INDEX root, std::vector<INDEX> &distances) {
    if (distances.size() != graph.nodes()){
        distances.clear();
        distances.resize(graph.nodes(), -1);
    }
    else{
        std::fill(distances.begin(), distances.end(), -1);
    }
    std::vector<bool> visitedNodes = std::vector<bool>(graph.nodes(), false);
    visitedNodes[root] = true;
    distances[root] = 0;
    std::deque<NodeId> nodes = std::deque<NodeId>();
    nodes.push_back(root);
    while (!nodes.empty()){
        NodeId currentNode = nodes.back();
        nodes.pop_back();
        for (NodeId i = 0; i < graph.degree(currentNode); ++i) {
            NodeId neighbor = graph.neighbor(currentNode, i);
            if (!visitedNodes[neighbor]){
                visitedNodes[neighbor] = true;
                distances[neighbor] = distances[currentNode] + 1;
                nodes.push_front(neighbor);
            }
        }
    }
}





template<typename T>
inline void GraphFunctions::DFS(const T &graph, T &tree, Nodes &nodeOrder, NodeId rootNodeId, const int seed) {
    std::mt19937_64 gen(seed);
    NodeId max = -1;
    if (rootNodeId == max) {
        rootNodeId = std::uniform_int_distribution<INDEX>(0, graph.nodes() - 1)(gen);
    }
    std::deque<std::pair<INDEX, INDEX>> swappedIds;
    std::vector<NodeId> neighborIds = std::vector<NodeId>(graph.maxDegree, 0);
    std::iota(neighborIds.begin(), neighborIds.end(), 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        tree.AddNodes(1);
    }
    tree.graphType = GraphType::TREE;
    std::unordered_set<INDEX> visitedNodes;
    std::vector<NodeId> CurrentNodes;
    CurrentNodes.emplace_back(rootNodeId);
    nodeOrder.emplace_back(rootNodeId);
    while (!CurrentNodes.empty()) {
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        if (visitedNodes.find(NextNodeId) == visitedNodes.end()) {
            visitedNodes.insert(NextNodeId);
            nodeOrder.emplace_back(NextNodeId);
            INDEX degree = graph.get_neighbors(NextNodeId).size();
            for (INDEX i = 0; i < degree; ++i) {
                //Get random neighbor
                INDEX idx = std::uniform_int_distribution<INDEX>(i, degree - 1)(gen);
                NodeId NeighborNodeId = graph.get_neighbors(NextNodeId)[neighborIds[idx]];
                std::swap(neighborIds[idx], neighborIds[i]);
                swappedIds.emplace_front(std::pair<int, int>{idx, i});
                if (visitedNodes.find(NeighborNodeId) == visitedNodes.end()) {
                    CurrentNodes.emplace_back(NeighborNodeId);
                }
            }
            tree.AddEdge(NextNodeId, CurrentNodes.back(), true);
            for (auto const &[a, b]: swappedIds) {
                std::swap(neighborIds[b], neighborIds[a]);
            }
            swappedIds.clear();
        }
    }
}

template<typename T>
inline void GraphFunctions::OptOrdering(const T &graph, Nodes &nodeOrder) {
    NodeId rootNodeId;
    NodeId max_degree = 0;
    for (int i = 0; i < graph.nodes(); ++i) {
        if (graph.degree(i) > max_degree){
            max_degree = graph.degree(i);
            rootNodeId = i;
        }
    }
    std::unordered_set<INDEX> visitedNodes = {rootNodeId};
    std::vector<NodeId> CurrentNodes;
    CurrentNodes.emplace_back(rootNodeId);
    nodeOrder.emplace_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        std::vector<std::pair<NodeId, INDEX>> node_degrees;
        for(auto n : graph.get_neighbors(NextNodeId)){
            node_degrees.emplace_back(n, graph.degree(n));
        }
        std::sort(node_degrees.begin(), node_degrees.end(), [&](const auto& a, const auto& b)
        {
            return a.second > b.second;
        });
        for (auto [n, d] : node_degrees) {
            if (visitedNodes.find(n) == visitedNodes.end()){
                CurrentNodes.emplace_back(n);
                nodeOrder.emplace_back(n);
                visitedNodes.insert(n);
            }
        }
    }
}

template<typename T>
inline void GraphFunctions::ReorderGraph(T &graph, const Nodes &nodeOrder) {

    std::unordered_map<NodeId, NodeId> nodeMap;
    int counter = 0;
    for (auto node : nodeOrder) {
        nodeMap.insert({node, counter});
        ++counter;
    }
    std::vector<Nodes> newNodes = std::vector<Nodes>(graph.nodes());
    counter = 0;
    for (auto node : nodeOrder) {
        const Nodes& neighbors = graph.get_neighbors(node);
        for (auto n : neighbors) {
            newNodes[counter].emplace_back(nodeMap[n]);
        }
        ++counter;
    }
    graph.set_graph(newNodes);

}

template<typename T1, typename T2>
inline bool GraphFunctions::CheckSpanningTree(const T1& graph, T2 &spanningTree) {
    if (spanningTree.nodes() != graph.nodes() || !GraphFunctions::CheckTree(spanningTree))
    {
        return false;
    }
    NodeId Counter = 0;
    for (auto const & Nodes : spanningTree.graph()) {
        for (auto Node : Nodes) {
            if (!graph.IsEdge(Counter, Node)) {
                return false;
            }
        }
        ++Counter;
    }
    return true;
}

template<typename T>
inline bool GraphFunctions::CheckTree(const T& graph) {
    for (INDEX i = 0; i < graph.nodes(); ++i) {
        if(graph.degree(i) == 0){
            return false;
        }
    }
    return graph.nodes() == graph.edges() + 1;
}






inline std::unordered_map<Label, Nodes> GraphFunctions::GetGraphLabelMap(Labels& nodeLabels) {
    std::unordered_map<Label, Nodes> graphLabelMap;
    for (NodeId i = 0; i < nodeLabels.size(); ++i) {
        if (graphLabelMap.find(nodeLabels[i]) == graphLabelMap.end()){
            graphLabelMap.insert(std::pair<Label, Nodes>(nodeLabels[i], std::vector<NodeId>{i}));
        }
        else{
            graphLabelMap.at(nodeLabels[i]).push_back(i);
        }
    }
    return graphLabelMap;
}

inline std::unordered_map<Label, INDEX> GraphFunctions::GetLabelFrequency(std::unordered_map<Label, Nodes> &labelMap) {
    std::unordered_map<Label, INDEX> labelFrequencyMap;
    for (auto const& [label, nodes] : labelMap){
        labelFrequencyMap.insert(std::pair<Label, INDEX>(label, nodes.size()));
    }
    return labelFrequencyMap;
}

inline std::pair<std::unordered_map<Label, Label>, std::unordered_map<Label, Label>> GraphFunctions::GetLabelMappingPair(Labels& nodeLabels) {
    std::unordered_map<Label, Label> labelMapping;
    std::unordered_map<Label, Label> reverseLabelMapping;
    std::set<Label> uniqueLabels = std::set<Label>{};
    for (auto label : nodeLabels) {
        uniqueLabels.insert(label);
    }
    for (auto label : uniqueLabels) {
        labelMapping.insert(std::pair<Label, Label>(label, labelMapping.size()));
        reverseLabelMapping.insert(std::pair<Label, Label>(labelMapping.size(), label));
    }
    return std::make_pair(labelMapping, reverseLabelMapping);
}

inline void GraphFunctions::GetNodesByLabel(Labels& nodeLabels, Label label, Nodes& nodesByLabel) {
    nodesByLabel.clear();
    for (Label i = 0; i < nodeLabels.size(); ++i) {
        if (nodeLabels[i] == label){
            nodesByLabel.push_back(i);
        }
    }
}

inline Labels GraphFunctions::GetLabelsByNodes(Nodes &nodeIds, const Labels &graphLabels) {
    Labels LabelsByNodes(nodeIds.size(), 0);
    for (INDEX i = 0; i < nodeIds.size(); ++i) {
        LabelsByNodes[i] = graphLabels[nodeIds[i]];
    }
    return LabelsByNodes;
}

inline std::unordered_set<Label> GraphFunctions::GetUniqueLabelsByNodes(Nodes &nodeIds, const Labels &GraphLabels) {
    std::unordered_set<Label> UniqueLabels = std::unordered_set<Label>{};
    Labels LabelsByNodes(nodeIds.size(), 0);
    for (auto nodeId : nodeIds) {
        UniqueLabels.insert(GraphLabels[nodeId]);
    }
    return UniqueLabels;
}



#endif //HOPS_GRAPHFUNCTIONS_H
