//
// Created by florian on 26.10.23.
//

#ifndef GOOGLE_TESTS_OUTERPLANARSUBGRAPHMITCHELL_H
#define GOOGLE_TESTS_OUTERPLANARSUBGRAPHMITCHELL_H

struct AlgorithmEdge{
public:
    AlgorithmEdge() = default;
    explicit AlgorithmEdge(const NodePair& nodePair, bool type = true);
    void setTreeEdge(bool b = true){ treeEdge = b;}
    [[nodiscard]] bool isTreeEdge() const{return treeEdge;};
    [[nodiscard]] bool is_valid() const{return _valid;};
    void set_valid(bool b){ _valid = b;};
    [[nodiscard]] int triangleCount() const{return _triangleCount;};
    [[nodiscard]] const NodePair& edge() const{return _nodePair;};
    void delete_edge(){_deleted = true; _valid = false; _triangleCount = 0;};
    [[nodiscard]] bool is_deleted() const {return _deleted;};
    void set_count(int count){_triangleCount = count;};
    bool original_edge() const;

private:
    NodePair _nodePair{};
    int _triangleCount = 0;
    bool _type = true;
    bool _valid = true;
    bool treeEdge = false;
    bool _deleted = false;
};

struct BiconnectedComponent{
    BiconnectedComponent()= default;
    //Nodes and degrees
    std::vector<NodeId> nodeIds;
    std::vector<NodeId> degree2Nodes;
    std::vector<int> nodeDegrees;

    //edges
    std::unordered_map<NodePair, AlgorithmEdge, hashNodePair> algorithmEdges;
    std::unordered_map<NodePair, bool, hashNodePair> edges;
    int componentId = -1;
};

class TriangulationForest{
public:
    enum class NodeType{
        EDGE,
        NODE,
    };
    struct TFNode{
        explicit TFNode(NodeType nodeType, AlgorithmEdge* correspondingEdge = nullptr, NodeId nodeId = -1)
                : nodeType(nodeType), correspondingEdge(correspondingEdge), nodeId(nodeId) {}

    public:
        std::vector<NodeId> children;
        NodeType nodeType;
        AlgorithmEdge* correspondingEdge;
        NodeId nodeId;
    };

    void AddNode(NodeType nodeType, AlgorithmEdge* edge = nullptr, NodeId nodeId = -1);
    void AddEdge(const NodePair& nodePair, NodeId nodeId);
    void AddEdge(NodeId nodeId, const NodePair& nodePair);
    int degree(const NodePair& nodePair){return (int) Node(nodePair)->children.size();};
    int degree(const AlgorithmEdge& algorithmEdge){ return degree(algorithmEdge.edge());};
    void clear(){_nodes.clear(); edgeMap.clear(); nodeMap.clear();}
    TFNode* rand_neighbor(TFNode* node, std::mt19937_64& gen);
    void GetValidPath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen);
    bool GetValidDeletePath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen);

    const TFNode* Node(NodeId nodeId) {return _node(nodeMap[nodeId]);};
    TFNode* GetNode(NodeId nodeId) {return _get_node(nodeMap[nodeId]);};
    const TFNode* Node(const NodePair& nodePair) {return _node(edgeMap[nodePair]);};
    TFNode* GetNode(const NodePair& nodePair) {return _get_node(edgeMap[nodePair]);};
    void print() const;

    void AddPath(std::vector<TFNode *> path);

private:
    std::vector<TFNode> _nodes;
    std::unordered_map<NodePair, int, hashNodePair> edgeMap;
    std::unordered_map<NodeId, int> nodeMap;
    TFNode* _get_node(NodeId nodeId){return &_nodes[nodeId];};
    const TFNode* _node(NodeId nodeId){return &_nodes[nodeId];};
};

class Triangle{
public:
    Triangle(const AlgorithmEdge& edgeA,const AlgorithmEdge& edgeB,const AlgorithmEdge& edgeBase) : a(edgeA), b(edgeB), base(edgeBase){};
    const AlgorithmEdge& a;
    const AlgorithmEdge& b;
    const AlgorithmEdge& base;
    void print() const{
        NodeId head = a.edge().first();
        if (a.edge().first() == base.edge().first() || a.edge().first() == base.edge().second()){
            head = a.edge().second();
        }
        std::cout << base.edge().first() <<  "/" << head << "\\" << base.edge().second() << " weights " << a.triangleCount() << "/" << base.triangleCount() << "\\" << b.triangleCount();
    }
};


class OuterPlanarSubgraphMitchell : public OuterplanarSubgraph {
public:
    explicit OuterPlanarSubgraphMitchell(const GraphStruct& graph);
    void generate(GraphStruct& subgraph, int seed, bool print) override;

    class PrintSteps {
    public:
        explicit PrintSteps(const TriangulationForest& triangulationForest) : _triangulationForest(triangulationForest){};
        NodePair currentDeleteEdge{};
        void printTriangulationStep(const Triangle& triangle) const;
        void printDeleteStep(std::vector<TriangulationForest::TFNode*>& path) const;
        void printRemovedNode(NodeId i) const;
        void printTriangulationForest() const;
        const TriangulationForest& _triangulationForest;
        void printEdge(const NodePair & pair) const;
    };
private:
    void getValidNeighbor(GraphStruct& graph, std::unordered_map<NodePair, AlgorithmEdge, hashNodePair> &algEdges,
                          int nodeId,
                          int neighborIdx, NodeId &neighborId);
    int nodeIdToComponentId(NodeId nodeId) const;
    void GetBiconnectedComponentsMitchell();
    void BiconnectedComponentSampling(GraphStruct& subgraph, BiconnectedComponent& biCom, std::mt19937_64& _gen);
    void Init();
    void UpdateDegrees(NodeId id);
    void DeleteRandomEdge(const PrintSteps& printing, std::mt19937_64& _gen);
    void CollapseTriangle(GraphStruct& component, AlgorithmEdge& edgeA, AlgorithmEdge& edgeB, NodePair& triangulationPair, const PrintSteps& printing, std::mt19937_64& _gen);
    void UpdateTriangulationForest(const NodePair &triangulationEdge, const NodePair &edgeA, const NodePair &edgeB, NodeId v);
    void DeleteEdges(std::vector<TriangulationForest::TFNode*>& path, const PrintSteps& printSteps, bool includingFirstEdge = false);


    std::vector<GraphStruct> biconnectedComponents;

    std::vector<std::list<int>> nodeToComponentAndId;
    std::vector<BiconnectedComponent> biComponents;
    BiconnectedComponent currentComponent;
    std::vector<NodePair> treeEdges;

    std::vector<NodePair> removableEdges;
    TriangulationForest triangulationForest;
    std::vector<TriangulationForest::TFNode*> _path1;
    std::vector<TriangulationForest::TFNode*> _path2;


    bool _print;
    int _seed;
};

OuterPlanarSubgraphMitchell::OuterPlanarSubgraphMitchell(const GraphStruct& graph)  : OuterplanarSubgraph(graph) {
    GetBiconnectedComponentsMitchell();
    Init();
}

AlgorithmEdge::AlgorithmEdge(const NodePair& nodePair, bool type)
        : _nodePair(nodePair), _type(type){}

bool AlgorithmEdge::original_edge() const {
    return _type;
}


void OuterPlanarSubgraphMitchell::getValidNeighbor(GraphStruct& graph, std::unordered_map<NodePair, AlgorithmEdge, hashNodePair>& algEdges, int nodeId, int neighborIdx, NodeId &neighborId) {
    neighborId = -1;
    int validIdx = -1;
    for (int i = 0; i < graph.degree(nodeId); ++i) {
        neighborId = graph.neighbor(nodeId, i);
        if(algEdges[NodePair(nodeId, neighborId)].is_valid()){
            ++validIdx;
        }
        if (validIdx == neighborIdx){
            return;
        }
    }
}

void OuterPlanarSubgraphMitchell::Init() {
    int componentCounter = 0;
    nodeToComponentAndId = std::vector<std::list<int>>(this->_graph.nodes(), std::list<int>());
    for (auto const & component : biconnectedComponents) {
        biComponents.emplace_back();
        BiconnectedComponent& biCom = biComponents.back();
        biCom.nodeIds = std::vector<NodeId>();
        biCom.nodeDegrees = std::vector<int>();
        biCom.algorithmEdges = std::unordered_map<NodePair, AlgorithmEdge, hashNodePair>();
        biCom.edges = std::unordered_map<NodePair, bool, hashNodePair>();
        biCom.degree2Nodes = std::vector<NodeId>();
        biCom.componentId = componentCounter;
        int componentNodeId = 0;
        for (auto nodeId = 0; nodeId < component.nodes(); ++nodeId) {
            nodeToComponentAndId[nodeId].push_back(componentNodeId);
            biCom.nodeIds.emplace_back(nodeId);
            biCom.nodeDegrees.emplace_back(component.degree(nodeId));
            if (component.degree(nodeId) == 2){
                biCom.degree2Nodes.emplace_back(nodeId);
            }
            ++componentNodeId;
        }

        // Iterate over all edges in the component and add them to the algorithm graph
        for (auto srcNode = 0; srcNode < component.nodes(); ++srcNode) {
            for (auto destNode = 0; destNode < component.degree(srcNode); ++destNode) {
                if (srcNode < destNode) {
                    NodePair nodePair = NodePair(srcNode, destNode);
                    biCom.algorithmEdges[nodePair] = AlgorithmEdge(nodePair);
                    biCom.edges[nodePair] = true;
                }
            }
        }
        ++componentCounter;
    }
}

void OuterPlanarSubgraphMitchell::GetBiconnectedComponentsMitchell() {
    std::vector<std::vector<NodeId>> components;
    GraphAlgorithms::GetBiconnectedComponents(_graph, components);
    for (auto component : components) {
        if (component.size() == 2){
            treeEdges.emplace_back(NodePair(component[0], component[1]));
        }
        else{
            this->biconnectedComponents.emplace_back(GraphStruct::SubGraph(_graph, component));
        }
    }

}

void OuterPlanarSubgraphMitchell::generate(GraphStruct& subgraph, int seed, bool p) {
    this->_seed = seed;
    this->_print = p;
    std::mt19937_64 _gen(seed);
    //Run spanning tree to mark non-deletable edges, TODO give spanning tree as parameter
    GraphStruct tree;
    std::vector<bool> visited;
    std::vector<INDEX> distances;
    BFSSpanningTree(_graph, tree, 0, visited, distances, false, seed);
    //GraphFunctions::bfsRandomSubtree(_graph, spanningTreeEdges, _gen);
    std::unordered_map<NodePair, bool, hashNodePair> edgesInSpanningTree;
    if (p){
        std::cout << "Spanning Tree: " << std::endl;
    }

    for (auto edge = tree.first_edge(); edge != tree.last_edge(); ++edge){
        edgesInSpanningTree[NodePair(*edge, false)] = true;
        if (p){
            NodePair(*edge, false).print();
        }
    }
    if (this->_print){
        std::cout << std::endl;
        std::cout  << std::endl << "\t" << "*************" << " Start Algorithm" << "******************" << std::endl;
    }
    for (BiconnectedComponent& biCom : biComponents) {
        for(auto & [nodePair, edgeInfo] : biCom.algorithmEdges){
            bool b = edgesInSpanningTree.find(nodePair) != edgesInSpanningTree.end();
            edgeInfo.setTreeEdge(b);
        }
        BiconnectedComponentSampling(subgraph, biCom, _gen);
        for (auto nodeId : biCom.nodeIds) {
            nodeToComponentAndId[nodeId].pop_front();
        }
    }
}

void OuterPlanarSubgraphMitchell::BiconnectedComponentSampling(GraphStruct& subgraph, BiconnectedComponent& biCom, std::mt19937_64& _gen) {
    GraphStruct& component = biconnectedComponents[biCom.componentId];
    INDEX nodeNum = component.nodes();
    currentComponent = biCom;

    //printing of the algorithm
    if (this->_print){
        std::cout  << std::endl << "\t" << "*************" << " New Component" << "******************" << std::endl;
    }

    removableEdges.clear();
    //Find all removable possibleEdges
    if (this->_print){
        std::cout << "\t" << "Removable Edges:" << std::endl << "\t";
    }
    for (auto const& [nodePair, edgeInfo]: currentComponent.algorithmEdges) {
        if (!edgeInfo.isTreeEdge()) {
            removableEdges.emplace_back(nodePair);
            if (this->_print){
                nodePair.print();
            }
        }
    }
    if (this->_print){
        std::cout << std::endl;
    }

    // iterate over all edges in the component
    for (auto edge = component.first_edge(); edge != component.last_edge(); ++edge) {
        triangulationForest.AddNode(TriangulationForest::NodeType::EDGE, &currentComponent.algorithmEdges[NodePair(*edge,false)]);
    }
    // iterate over all nodes in the component
    for (unsigned int node : component) {
        triangulationForest.AddNode(TriangulationForest::NodeType::NODE, nullptr, node);
    }
    auto printing = PrintSteps(triangulationForest);
    while (nodeNum > 0){
        std::vector<NodeId>& degree2Nodes = currentComponent.degree2Nodes;
        if (this->_print){
            std::cout << "\t" << "Next Step: " << std::endl;
            std::cout << "\t\t Node Degrees: ";
            for (int i = 0; i < currentComponent.nodeIds.size(); ++i) {
                std::cout << "(" << currentComponent.nodeIds[i] << ", " << currentComponent.nodeDegrees[i] << " deg)";
            }
            std::cout << std::endl;
            std::cout << "\t\t Degree 2 Nodes: ";
            for (auto nodeId : degree2Nodes) {
                std::cout << nodeId << ", ";
            }
            std::cout << std::endl;
        }
        if (!degree2Nodes.empty()){
            int rnd = std::uniform_int_distribution<int>(0, (int) degree2Nodes.size() - 1)(_gen);
            NodeId v = degree2Nodes[rnd];
            std::swap(degree2Nodes[rnd], degree2Nodes.back());
            degree2Nodes.pop_back();
            if (currentComponent.nodeDegrees[nodeIdToComponentId(v)] == 2){
                currentComponent.nodeDegrees[nodeIdToComponentId(v)] = 0;
                if (this->_print){
                    std::cout << "\t\t\t" << "Triangulate with get_node " << v << std::endl;
                }
                --nodeNum;
                NodeId x, y;
                getValidNeighbor(component, currentComponent.algorithmEdges, v, 0, x);
                getValidNeighbor(component, currentComponent.algorithmEdges, v, 1, y);
                NodePair edgeA = NodePair(v, x);
                NodePair edgeB = NodePair(v, y);
                AlgorithmEdge& algorithmEdgeA = currentComponent.algorithmEdges[edgeA];
                AlgorithmEdge& algorithmEdgeB = currentComponent.algorithmEdges[edgeB];
                NodePair triangulationEdge = NodePair(x, y);
                CollapseTriangle(component, algorithmEdgeA, algorithmEdgeB, triangulationEdge, printing, _gen);
                UpdateTriangulationForest(triangulationEdge, edgeA, edgeB, v);
            }
            else{
                --nodeNum;
                if (currentComponent.nodeDegrees[nodeIdToComponentId(v)] == 1){
                    NodeId x;
                    getValidNeighbor(component, currentComponent.algorithmEdges, v, 0, x);
                    --currentComponent.nodeDegrees[nodeIdToComponentId(x)];
                    if (currentComponent.nodeDegrees[nodeIdToComponentId(x)] == 0){
                        --nodeNum;
                    }
                }
                currentComponent.nodeDegrees[nodeIdToComponentId(v)] = 0;
                if (_print){
                    printing.printRemovedNode(v);
                }
            }

        }
        else {
            DeleteRandomEdge(printing, _gen);
        }
        if (this->_print) {
            printing.printTriangulationForest();
        }
    }
    //Add all undeleted edges to the _subgraph
    for (auto const& [edge, valid] : currentComponent.edges) {
        if (valid) {
           subgraph.AddEdge(edge.first(), edge.second());
        }
    }
}

void OuterPlanarSubgraphMitchell::DeleteRandomEdge(const PrintSteps& printing, std::mt19937_64& _gen) {
    while(!removableEdges.empty()){
        //Get random edge from all the edges which can be removed
        int idx = std::uniform_int_distribution<int>(0, (int) removableEdges.size() - 1)(_gen);
        NodePair nodePair = removableEdges[idx];
        std::swap(removableEdges[idx], removableEdges.back());
        removableEdges.pop_back();

        //Get corresponding algorithm edge and check if it is valid
        AlgorithmEdge& algorithmEdge = currentComponent.algorithmEdges[nodePair];
        if (algorithmEdge.is_valid()) {
            auto src = nodePair.first(), dst = nodePair.second();
            auto componentSrc = nodeIdToComponentId(src), componentDst = nodeIdToComponentId(dst);
            int degree = triangulationForest.degree(algorithmEdge);
            if (degree == algorithmEdge.triangleCount()){
                bool validPath = triangulationForest.GetValidDeletePath(nodePair, _path1, _gen);
                if (validPath) {
                    if (degree == 2) {
                        if (triangulationForest.GetValidDeletePath(nodePair, _path2, _gen) && _path2.size() > 1) {
                            if (_print){
                                printing.printEdge(nodePair);
                            }
                            DeleteEdges(_path1, printing, true);
                            DeleteEdges(_path2, printing, true);
                            --currentComponent.nodeDegrees[componentSrc];
                            --currentComponent.nodeDegrees[componentDst];
                            if (currentComponent.nodeDegrees[componentSrc] == 2) {
                                currentComponent.degree2Nodes.emplace_back(src);
                            }
                            if (currentComponent.nodeDegrees[componentDst] == 2) {
                                currentComponent.degree2Nodes.emplace_back(dst);
                            }
                            break;
                        }
                        else{
                            triangulationForest.AddPath(_path1);
                        }
                    } else {
                        if (_print){
                            printing.printEdge(nodePair);
                        }
                        DeleteEdges(_path1, printing, true);
                        --currentComponent.nodeDegrees[componentSrc];
                        --currentComponent.nodeDegrees[componentDst];
                        if (currentComponent.nodeDegrees[componentSrc] == 2) {
                            currentComponent.degree2Nodes.emplace_back(src);
                        }
                        if (currentComponent.nodeDegrees[componentDst] == 2) {
                            currentComponent.degree2Nodes.emplace_back(dst);
                        }
                        break;
                    }
                }
            }
        }
    }
}

void OuterPlanarSubgraphMitchell::UpdateDegrees(NodeId id) {
    int componentId = nodeIdToComponentId(id);
    --currentComponent.nodeDegrees[componentId];
    if (currentComponent.nodeDegrees[componentId] == 2){
        currentComponent.degree2Nodes.emplace_back(id);
    }
}

void OuterPlanarSubgraphMitchell::CollapseTriangle(GraphStruct& component, AlgorithmEdge& edgeA, AlgorithmEdge& edgeB, NodePair& triangulationPair, const PrintSteps& printing, std::mt19937_64& _gen) {
    AlgorithmEdge* triangulationEdge = nullptr;
    edgeA.set_valid(false);
    edgeB.set_valid(false);
    if (currentComponent.algorithmEdges.find(triangulationPair) != currentComponent.algorithmEdges.end()){
        triangulationEdge = &currentComponent.algorithmEdges[triangulationPair];
    }

    if (edgeA.triangleCount() > 1){
        triangulationForest.GetValidPath(edgeA.edge(), _path1, _gen);
        DeleteEdges(_path1, printing);
        edgeA.set_count(edgeA.triangleCount() - 1);
    }
    if (edgeB.triangleCount() > 1){
        triangulationForest.GetValidPath(edgeB.edge(), _path1, _gen);
        DeleteEdges(_path1, printing);
        edgeB.set_count(edgeB.triangleCount() - 1);
    }
    if (triangulationEdge != nullptr && triangulationEdge->triangleCount() > 1){
        triangulationForest.GetValidPath(triangulationEdge->edge(), _path1, _gen);
        DeleteEdges(_path1, printing);
        triangulationEdge->set_count(triangulationEdge->triangleCount() - 1);
    }

    //Triangulated baseEdge is not contained in the graph
    if (triangulationEdge == nullptr || !triangulationEdge->is_valid()){
        if (triangulationEdge == nullptr){
            component.AddEdge(triangulationPair.first(), triangulationPair.second());
            currentComponent.algorithmEdges[triangulationPair] = AlgorithmEdge(triangulationPair, false);
            triangulationEdge = &currentComponent.algorithmEdges[triangulationPair];
            triangulationForest.AddNode(TriangulationForest::NodeType::EDGE, triangulationEdge);
        }
        triangulationEdge->set_valid(true);
        triangulationEdge->setTreeEdge(false);
        triangulationEdge->set_count(1);
        removableEdges.emplace_back(triangulationPair);
    }
        //Triangulation edge is in the graph
    else{
        triangulationEdge->set_count(triangulationEdge->triangleCount() + 1);
        UpdateDegrees(triangulationPair.first());
        UpdateDegrees(triangulationPair.second());
    }

    if (this->_print) {
        printing.printTriangulationStep(Triangle(edgeA, edgeB, *triangulationEdge));
    }
}

void OuterPlanarSubgraphMitchell::UpdateTriangulationForest(const NodePair &triangulationEdge, const NodePair &edgeA,
                                                            const NodePair &edgeB, NodeId v) {

    //Add collapsed edges to collapse tree and connect to base edge
    triangulationForest.AddEdge(triangulationEdge, v);
    triangulationForest.AddEdge(v, edgeA);
    triangulationForest.AddEdge(v, edgeB);

}

void OuterPlanarSubgraphMitchell::DeleteEdges(std::vector<TriangulationForest::TFNode*> &path, const PrintSteps& printSteps, bool includingFirstEdge) {
    for (auto const node : path) {
        if (!includingFirstEdge){
            includingFirstEdge = true;
        }
        else {
            if (node->nodeType == TriangulationForest::NodeType::EDGE) {
                node->correspondingEdge->delete_edge();
                auto const it = currentComponent.edges.find(node->correspondingEdge->edge());
                if (it != currentComponent.edges.end()){
                    currentComponent.edges[node->correspondingEdge->edge()] = false;
                }
            }
        }
    }
    if (this->_print) {
        printSteps.printDeleteStep(path);
    }
}


void TriangulationForest::AddNode(TriangulationForest::NodeType nodeType, AlgorithmEdge* edge, NodeId nodeId) {
    _nodes.emplace_back(TFNode(nodeType, edge, nodeId));
    if (nodeType == NodeType::EDGE){
        edgeMap[edge->edge()] = (int) _nodes.size() - 1;
    }
    else{
        nodeMap[nodeId] = (int) _nodes.size() - 1;
    }
}

void TriangulationForest::AddEdge(NodeId nodeId, const NodePair& nodePair) {
    TFNode* srcNode = GetNode(nodeId);
    srcNode->children.emplace_back(edgeMap[nodePair]);
}

void TriangulationForest::AddEdge(const NodePair &nodePair, NodeId nodeId) {
    TFNode* srcNode = GetNode(nodePair);
    srcNode->children.emplace_back(nodeMap[nodeId]);
}

TriangulationForest::TFNode* TriangulationForest::rand_neighbor(TFNode* node, std::mt19937_64 &gen) {
    TFNode* neighbor = nullptr;
    if (!node->children.empty()){
        int rndIdx = std::uniform_int_distribution<int>(0, (int) node->children.size() - 1)(gen);
        neighbor = _get_node(node->children[rndIdx]);
        node->children.erase(node->children.begin() + rndIdx);
    }
    return neighbor;
}

void TriangulationForest::GetValidPath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen) {
    path.clear();
    TFNode* pathRoot = GetNode(nodePair);
    path.emplace_back(GetNode(nodePair));
    bool pathUp = false;
    while (!path.empty()){
        TFNode* neighbor = rand_neighbor(path.back(), gen);
        if (neighbor == nullptr){
            if (!pathUp) {
                break;
            }
            else{
                path.pop_back();
            }
        }
        else{
            if (neighbor->correspondingEdge == nullptr || !neighbor->correspondingEdge->isTreeEdge()){
                path.emplace_back(neighbor);
                pathUp = false;
            }
            else{
                pathUp = true;
            }
        }
    }
    if (pathRoot->children.empty()){
        pathRoot->correspondingEdge->setTreeEdge();
    }
}

bool TriangulationForest::GetValidDeletePath(const NodePair &nodePair, std::vector<TFNode*>& path, std::mt19937_64& gen) {
    GetValidPath(nodePair, path, gen);
    return (int) path.size() > 0 && path.back()->correspondingEdge->original_edge();
}

void TriangulationForest::print() const{
    for (auto const & node : _nodes) {
        for (auto const & children : node.children) {
            if (node.nodeType == NodeType::NODE){
                std::cout << node.nodeId;
            }
            else{
                std::cout << "(";
                node.correspondingEdge->edge().print();
                std::cout << ")";
            }
            std::cout << " - ";
            const TFNode& childNode = _nodes[children];
            if (childNode.nodeType == NodeType::NODE){
                std::cout << childNode.nodeId << " ";
            }
            else{
                std::cout << "(";
                childNode.correspondingEdge->edge().print();
                std::cout << ") ";
            }
        }
    }
    std::cout << std::endl;
}

void TriangulationForest::AddPath(std::vector<TFNode *> path) {
    for (int i = 0; i < path.size() - 1; ++i) {
        if (path[i+1]->nodeType == NodeType::NODE){
            path[i]->children.emplace_back(nodeMap[path[i+1]->nodeId]);
        }
        else{
            path[i]->children.emplace_back(edgeMap[path[i+1]->correspondingEdge->edge()]);
        }
    }
}


void OuterPlanarSubgraphMitchell::PrintSteps::printTriangulationStep(const Triangle& triangle) const{
    std::cout << "\t\t\t\t" << "Triangulate ";
    triangle.print();
    std::cout << std::endl;
}

void OuterPlanarSubgraphMitchell::PrintSteps::printDeleteStep(std::vector<TriangulationForest::TFNode*>& path) const {
    std::cout << "\t\t\t\t" << "Delete Edges ";
    for (auto const & node : path) {
        if (node->nodeType == TriangulationForest::NodeType::EDGE){
            std::cout << " ";
            node->correspondingEdge->edge().print();
            std::cout << " ";
        }
    }
    std::cout << std::endl;
}

void OuterPlanarSubgraphMitchell::PrintSteps::printRemovedNode(NodeId i) const {
    std::cout << "\t\t\t" << "Removed Node " << i;
    std::cout << std::endl;
}

void OuterPlanarSubgraphMitchell::PrintSteps::printTriangulationForest() const {
    std::cout << "\t\t\t\t" << "New Triangulation Forest ";
    _triangulationForest.print();
}

void OuterPlanarSubgraphMitchell::PrintSteps::printEdge(const NodePair & pair) const {
    std::cout << "\t\t\t\t" << "Delete Base Edge ";
    pair.print();
}

int OuterPlanarSubgraphMitchell::nodeIdToComponentId(NodeId nodeId) const {
    return nodeToComponentAndId[nodeId].front();
}



#endif //GOOGLE_TESTS_OUTERPLANARSUBGRAPHMITCHELL_H
