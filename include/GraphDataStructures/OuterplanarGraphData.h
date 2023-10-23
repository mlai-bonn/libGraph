//
// Created by anonymous on 07.11.21.
//

#ifndef CLOSURES_OUTERPLANARGRAPHDATA_H
#define CLOSURES_OUTERPLANARGRAPHDATA_H


#include "Algorithms/Graph/GraphAlgorithms.h"

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

class OuterplanarGraphData : public GraphStruct {
public:
    OuterplanarGraphData() : GraphStruct(){};
    explicit OuterplanarGraphData(const std::string & path) : GraphStruct(path){
        nodeToComponents = std::vector<std::vector<int>>(this->nodes());
        maxCompSize = std::vector<int>(this->nodes(), 0);
        set();
    };
    OuterplanarGraphData(GraphStruct& pGraph, int size) : GraphStruct(size, {}){
        nodeToComponents = std::vector<std::vector<int>>(pGraph.nodes());
        maxCompSize = std::vector<int>(pGraph.nodes(), 0);
        set();
    };

    std::vector<OuterplanarComponent> Components;
    void set();
    GraphStruct& get_bbTree() {return bbTree.tree;};
    NodeOrComponent& get_bbNodeOrComponent(NodeId nodeId){return bbTree.BBNodesToNodeOrComponent[nodeId];};
    void get_components(NodeId nodeId, std::vector<NodeOrComponent>& nodeOrComponent);
    int get_component_num(NodeId nodeId){return (int) bbTree.nodeToBBNodes[nodeId].size();};
    int max_component_size(NodeId nodeId){return maxCompSize[nodeId];};
    void get_bb_tree_ids(NodeId nodeId, std::set<NodeId>& ids);

private:
    std::vector<std::vector<int>> nodeToComponents;
    std::vector<int> maxCompSize;
    BBTree bbTree;

    void init_outerplanar();
};

void OuterplanarGraphData::set() {
    if (nodeToComponents.empty()){
        this->init_outerplanar();
    }
    bbTree.tree.set_type(GraphType::TREE);
    std::vector<std::vector<NodeId>> components;
    GetBiconnectedComponents(static_cast<GraphStruct>(*this), components);
    //StaticFunctions::printComponents(components);
    for (int i = 0; i < components.size(); ++i) {
        auto & currentComponent = components[i];
        if (currentComponent.size() > 2) {
            Components.emplace_back();
            Components.back().component = SubGraph(static_cast<GraphStruct>(*this), currentComponent);
            Components.back().component.set_type(GraphType::OUTERPLANAR);
            Components.back().getInnerFaces();
        }
        for (int j = 0; j < currentComponent.size(); ++j) {
            NodeId currentNodeId = currentComponent[j];
            nodeToComponents[currentNodeId].emplace_back(i);
            maxCompSize[currentNodeId] = std::max(maxCompSize[currentNodeId], (int) currentComponent.size());
            if (currentComponent.size() > 2) {
                Components.back().NodeIdToComponentNodeId[currentNodeId] = j;
                Components.back().ComponentNodeIdToNodeId[j] = currentNodeId;
            }
        }
    }
    for (int i = 0; i < nodeToComponents.size(); ++i) {
        auto & comps = nodeToComponents[i];
        if (comps.size() > 1 || max_component_size(i) == 2){
            NodeId currentNode = bbTree.tree.add_node();
            bbTree.nodeToBBNodes[i] = {currentNode};
            bbTree.BBNodesToNodeOrComponent[currentNode] = {i, -1};
        }
        else{
            bbTree.nodeToBBNodes[i] = {};
        }
    }
    int c = 0;
    for (auto & component : components) {
        auto & currentComponent = component;
        if (currentComponent.size() > 2) {
            NodeId component_node = bbTree.tree.add_node();
            for (unsigned int j : currentComponent) {
                bbTree.nodeToBBNodes[j].emplace_back(component_node);
                bbTree.BBNodesToNodeOrComponent[component_node] = {-1, c};
                NodeId currentNodeId = j;
                if (nodeToComponents[currentNodeId].size() > 1){
                    bbTree.tree.add_edge(bbTree.nodeToBBNodes[currentNodeId][0], component_node);
                }
            }
        }
        else{
            --c;
            bbTree.tree.add_edge(bbTree.nodeToBBNodes[component[0]][0], bbTree.nodeToBBNodes[component[1]][0]);
        }
        ++c;
    }
    //bbTree.tree.print();
    //bbTree.tree.init();
}

void OuterplanarGraphData::get_components(NodeId nodeId, std::vector<NodeOrComponent> &nodeOrComponent) {
    nodeOrComponent.clear();
    for(auto const & bbnode : bbTree.nodeToBBNodes[nodeId]){
        nodeOrComponent.emplace_back(bbTree.BBNodesToNodeOrComponent[bbnode]);
    }
}

void OuterplanarGraphData::get_bb_tree_ids(NodeId nodeId, std::set<NodeId>& ids) {
    ids.insert(bbTree.nodeToBBNodes[nodeId].begin(), bbTree.nodeToBBNodes[nodeId].end());
}

void OuterplanarGraphData::init_outerplanar() {
    nodeToComponents = std::vector<std::vector<int>>(nodes());
    maxCompSize = std::vector<int>(nodes(), 0);
}

void OuterplanarComponent::getInnerFaces() {

}

BBTree::~BBTree() {
}



#endif //CLOSURES_OUTERPLANARGRAPHDATA_H
