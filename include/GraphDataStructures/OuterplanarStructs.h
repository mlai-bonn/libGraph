//
// Created by florian on 24.10.23.
//

#ifndef GOOGLE_TESTS_OUTERPLANARSTRUCTS_H
#define GOOGLE_TESTS_OUTERPLANARSTRUCTS_H

#include <map>

struct OuterplanarComponent{
    std::vector<GraphStruct> faces;
    std::vector<std::vector<int>> nodeToFaces;
    GraphStruct component = GraphStruct();
    std::map<NodeId, NodeId> NodeIdToComponentNodeId;
    std::map<NodeId, NodeId> ComponentNodeIdToNodeId;
    void getInnerFaces();
    NodeId componentId(NodeId nodeId){return NodeIdToComponentNodeId[nodeId];};
    NodeId nodeId(NodeId componentNodeId){return ComponentNodeIdToNodeId[componentNodeId];};
};



struct NodeOrComponent{
    int _nodeId = -1;
    int _componentId = -1;

    bool is_node(int& id) const{id = _nodeId; return _nodeId != -1;};
    bool is_component(int& id) const{id = _componentId; return _componentId != -1;};
};

struct BBTree{
    ~BBTree();
    GraphStruct tree = GraphStruct();
    std::unordered_map<NodeId, std::vector<NodeId>> nodeToBBNodes;
    std::unordered_map<NodeId, NodeOrComponent> BBNodesToNodeOrComponent;
};


#endif //GOOGLE_TESTS_OUTERPLANARSTRUCTS_H
