//
// Created by Florian on 10.03.2021.
//

#ifndef HOPS_ROOTEDPATTERN_H
#define HOPS_ROOTEDPATTERN_H

#include <vector>
#include <map>

class RootedPattern {
public:
    RootedPattern() = default;

    explicit RootedPattern(DGraphStruct& tree, NodeId RootNode = -1, LABEL_TYPE labelType = LABEL_TYPE::UNLABELED);
    void GetBFSOrder(NodeId RootNodeId);
    NodeId GetBFSOrderIndexByNodeId(NodeId Id) const;

    DGraphStruct tree;
    NodeId RootNode = 0;
    Nodes BFSOrderIndex;
    size_t NumberOfPartitions = 0;
    std::vector<std::pair<NodeId, Nodes>> BFSOrder = std::vector<std::pair<NodeId, Nodes>>();
    std::vector<std::map<Label, Nodes>> BFSOrderLabelMap = std::vector<std::map<Label, Nodes>>();
};

inline RootedPattern::RootedPattern(DGraphStruct& tree, NodeId RootNode, LABEL_TYPE labelType) : tree(tree), RootNode(RootNode) {
    this->BFSOrderIndex = Nodes(this->tree.nodes(), std::numeric_limits<INDEX>::max());
    if (this->RootNode != std::numeric_limits<INDEX>::max()) {
        this->GetBFSOrder(this->RootNode);
        if (labelType != LABEL_TYPE::UNLABELED && !this->tree.labels().empty()) {
            for (std::pair<NodeId, Nodes> &ParentChildrenPair : this->BFSOrder) {
                std::unordered_set < Label > UniqueLabels = GraphFunctions::GetUniqueLabelsByNodes(ParentChildrenPair.second,
                                                                                                   this->tree.labels());
                this->BFSOrderLabelMap.emplace_back();
                for (const Label &label : UniqueLabels) {
                    this->BFSOrderLabelMap.back().insert(std::pair<Label, Nodes>(label, Nodes()));
                }
                Labels ChildrenLabels = GraphFunctions::GetLabelsByNodes(ParentChildrenPair.second, this->tree.labels());
                for (size_t i = 0; i < ChildrenLabels.size(); ++i) {
                    this->BFSOrderLabelMap.back()[ChildrenLabels[i]].push_back(ParentChildrenPair.second[i]);
                }
            }
        }
    }
}



inline void RootedPattern::GetBFSOrder(NodeId RootNodeId) {
    BFSOrder.clear();
    std::vector<bool> VisitedNodes(tree.nodes(), false);
    std::vector<NodeId> CurrentNodes = std::vector<NodeId>{RootNodeId};
    VisitedNodes[RootNodeId] = true;
    NodeId NodeCounter = 0;
    BFSOrderIndex[RootNodeId] = 0;
    NodeId Current_Node;
    while (!CurrentNodes.empty()) {
        Current_Node = CurrentNodes.back();
        CurrentNodes.pop_back();
        NodeId Degree = tree.out_degree(Current_Node);
        if (Degree > 0) {
            BFSOrder.emplace_back(Current_Node, Nodes());
            ++NumberOfPartitions;
            for (int Id = 0; Id < Degree; ++Id) {
                NodeId Child_Node = tree.neighbor(Current_Node, Id);
                if (!VisitedNodes[Child_Node]) {
                    CurrentNodes.insert(CurrentNodes.begin(), Child_Node);
                    ++NodeCounter;
                    BFSOrderIndex[Child_Node] = NodeCounter;
                    BFSOrder.back().second.push_back(Child_Node);
                    VisitedNodes[Child_Node] = true;
                }
            }
        }
    }
}

inline NodeId RootedPattern::GetBFSOrderIndexByNodeId(NodeId Id) const {
    return BFSOrderIndex[Id];
}

#endif //HOPS_ROOTEDPATTERN_H
