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


struct GraphStruct;
class DGraphStruct;
class DDataGraph;

class GraphFunctions {
public:
    static std::unordered_map<Label, INDEX> GetLabelFrequency(std::unordered_map<Label, Nodes>& labelMap);
    static std::unordered_map<Label, Nodes> GetGraphLabelMap(Labels& nodeLabels);
    static std::pair<std::unordered_map<Label, Label>, std::unordered_map<Label, Label>> GetLabelMappingPair(Labels& nodeLabels);
    static void GetNodesByLabel(Labels& nodeLabels, Label label, Nodes& nodesByLabel);
    static Labels GetLabelsByNodes(Nodes& nodeIds, const Labels& graphLabels);
    static std::unordered_set<Label> GetUniqueLabelsByNodes(Nodes& nodeIds, const Labels& GraphLabels);
};

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
