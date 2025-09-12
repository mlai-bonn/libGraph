//
// Created by anonymous on 10.06.2021.
//

#ifndef CLOSURES_OUTERPLANARSUBGRAPHDFS_H
#define CLOSURES_OUTERPLANARSUBGRAPHDFS_H


#include "OuterplanarSubgraph.h"

enum class Side {
    LEFT,
    RIGHT,
};

class NodeReachability{
public:
    NodeReachability(): _left(true), _right(true){};
    [[nodiscard]] bool left() const {return _left;};
    [[nodiscard]] bool right() const {return _right;};
    void set(bool left, bool right){ _left = left; _right = right;};
    void close_left(){ _left = false;};
    void open_left(){_left = true;};
    void close_right(){ _right = false;};
    void open_right(){_right = true;};
    void reset(){ _right = true, _left = true;};
    [[nodiscard]] int count() const{
        int count = 0;
        if (_left){++count;};
        if(_right){++count;};
        return count;
    }
private:
    bool _left = true;
    bool _right = true;

};

class NodeStruct{
public:
    NodeStruct() = default;
    explicit NodeStruct(int id) : node_id(id){};
    NodeId node_id = -1;
    bool visited = false;
    int dfs_depth = -1;
    NodeReachability nodeReachability = NodeReachability();
    NodeStruct* dfs_parent = nullptr;
    int last_R = 0;
    int last_L = 0;

    void reset(){
        dfs_parent = nullptr;
        visited = false;
        dfs_depth = -1;
        nodeReachability.reset();
        last_L = 0;
        last_R = 0;
        left_spanning = false;
        right_spanning = false;
    };

    bool operator <(const NodeStruct &b) const{
        return this->dfs_depth < b.dfs_depth;
    }
    bool operator >(const NodeStruct &b) const{
        return this->dfs_depth > b.dfs_depth;
    }
    bool operator ==(const NodeStruct &b) const{
        return this->dfs_depth == b.dfs_depth;
    }
    bool operator >(int i) const{
        return this->dfs_depth > i;
    }
    bool operator <(int i) const{
        return this->dfs_depth < i;
    }
    bool operator >=(int i) const{
        return this->dfs_depth >= i;
    }
    bool operator <=(int i) const{
        return this->dfs_depth <= i;
    }
    bool operator ==(int i) const{
        return this->dfs_depth == i;
    }

    NodeStruct *last_spanned_root = nullptr;
    bool left_spanning = false;
    bool right_spanning = false;
};

class DirectedEdgeStruct{
public:
    DirectedEdgeStruct(NodeId src, NodeId dst) : src(src), dst(dst){};
    DirectedEdgeStruct(const NodeStruct& src, const NodeStruct& dst) : src(src.node_id), dst(dst.node_id){};
    NodeId src;
    NodeId dst;
};


class OuterplanarSubgraphDFS : public OuterplanarSubgraph {
public:
    explicit OuterplanarSubgraphDFS(const GraphStruct& graph);

    //Copy constructor
    OuterplanarSubgraphDFS(const OuterplanarSubgraphDFS& other);

    void generate(GraphStruct& graphStruct, int seed, bool p) override;
    void get_next_node(GraphStruct& subgraph, std::mt19937_64& gen);
    static bool CheckLeft(const NodeStruct& src, const NodeStruct& dst) ;
    static bool CheckRight(const NodeStruct& src, const NodeStruct& dst) ;
    [[nodiscard]] static bool crossesLeftDiagonal(const NodeStruct& dst) ;
    [[nodiscard]] static bool crossesRightDiagonal(const NodeStruct& dst) ;
    [[nodiscard]] static bool enclosingPointByLeftEdge(const NodeStruct &src,const NodeStruct &dst) ;
    [[nodiscard]] static bool enclosingPointByRightEdge(const NodeStruct &src,const NodeStruct &dst) ;

    void GetCotreeEdges(NodeStruct& currentNode, std::mt19937_64& gen);

    void GetValidEdges(NodeStruct& currentNode, std::vector<DirectedEdgeStruct>& edges, Side side);
    void AddEdges(const std::vector<DirectedEdgeStruct>& edges, Side side);
    void AddGraphEdges(GraphStruct& subgraph, const std::vector<DirectedEdgeStruct>& edges, Side side);

    void reset();

    std::vector<NodeId> reachableNodes;
private:
    std::vector<NodeStruct> graphNodes;
    NodeId dfs_root_node{};
    std::vector<DirectedEdgeStruct> possibleEdges;
    std::vector<NodeId> neighborsToVisit;
    std::vector<NodeId> dfs_stack;
    std::vector<DirectedEdgeStruct> left_edges;
    std::vector<DirectedEdgeStruct> right_edges;
    //Used to take neighbor uniformly at random
    std::vector<int> neighborIds;
    std::vector<int> swaps;
    bool print = false;
    bool dfs_tree = false;

    void UpdateNodeParameters(NodeStruct &currentNode);
};

OuterplanarSubgraphDFS::OuterplanarSubgraphDFS(const GraphStruct & graph) : OuterplanarSubgraph(graph) {
    INDEX maxDegree = 0;
    graphNodes = std::vector<NodeStruct>(graph.nodes(), NodeStruct());
    for (int i = 0; i < graph.nodes(); ++i) {
        graphNodes[i] = NodeStruct(i);
        maxDegree = std::max(maxDegree, graph.degree(i));
    }
    neighborIds = std::vector<int>(maxDegree);
    std::iota(neighborIds.begin(), neighborIds.end(), 0);
}

void OuterplanarSubgraphDFS::generate(GraphStruct& graphStruct, int seed, bool p) {
    graphStruct.Reset(_graph.nodes());
    graphStruct.SetType(GraphType::OUTERPLANAR);
    std::mt19937_64 _gen(seed);
    this->print = p;
    this->reset();
    dfs_root_node = std::uniform_int_distribution<INDEX>(0, this->_graph.nodes() - 1)(_gen);
    dfs_stack.emplace_back(dfs_root_node);
    NodeStruct* path_root = &this->graphNodes[dfs_root_node];
    path_root->last_spanned_root = path_root;
    path_root->dfs_parent = path_root;
    path_root->dfs_depth = 0;

    if (p) {
        std::cout << std::endl;
    }
    while (!this->dfs_stack.empty()) {
        get_next_node(graphStruct, _gen);
    }
    if (p) {
        std::cout << std::endl;
    }
}

void OuterplanarSubgraphDFS::get_next_node(GraphStruct& subgraph, std::mt19937_64 & gen) {
    //Current _node in the dfs
    NodeId currentNodeId = this->dfs_stack.back();
    dfs_stack.pop_back();
    NodeStruct &currentNode = this->graphNodes[currentNodeId];
    //If current get_node is not visited go further down in the dfs and update parameters
    if (!currentNode.visited) {
        //Check if get_node is not the root node
        if (currentNodeId != dfs_root_node) {
            //Add tree edge to the _graph
            subgraph.AddEdge(currentNodeId, currentNode.dfs_parent->node_id);
            UpdateNodeParameters(currentNode);
        }

        //Set get_node to visited
        currentNode.visited = true;
        reachableNodes.emplace_back(currentNodeId);

        //Traverse all the neighbors of currentNode and get cotree edges
        GetCotreeEdges(currentNode, gen);
        if (!neighborsToVisit.empty()) {
            GetValidEdges(currentNode, left_edges, Side::LEFT);
            GetValidEdges(currentNode, right_edges, Side::RIGHT);
            if (left_edges.size() >= right_edges.size()) {
                AddEdges(left_edges, Side::LEFT);
                AddGraphEdges(subgraph, left_edges, Side::LEFT);
            } else {
                AddEdges(right_edges, Side::RIGHT);
                AddGraphEdges(subgraph, right_edges, Side::RIGHT);
            }
        }
    }
}


void OuterplanarSubgraphDFS::GetCotreeEdges(NodeStruct& currentNode, std::mt19937_64& gen) {
    //Traverse all neighbors of get_node
    //auto node = this->_graph->GetNI(currentNode.node_id);
    INDEX degree = this->_graph.degree(currentNode.node_id);
    this->neighborsToVisit.clear();
    INDEX randIdx;
    NodeId neighbor;
    NodeStruct* neighborNode;
    for (int i = 0; i < degree; ++i) {
        randIdx = std::uniform_int_distribution<INDEX>(i, degree - 1)(gen);
        neighbor = this->_graph.neighbor(currentNode.node_id, neighborIds[randIdx]);
        swaps.emplace_back(randIdx);
        std::swap(neighborIds[i], neighborIds[randIdx]);
        if (!currentNode.dfs_parent || neighbor != currentNode.dfs_parent->node_id) {
            if (this->graphNodes[neighbor].visited) {
                //Edge currentNodeId, visitedNeighbor can be some diagonal
                this->neighborsToVisit.emplace_back(neighbor);
            } else {
                //Add node to dfs_stack and update parent get_node
                neighborNode = &this->graphNodes[neighbor];
                neighborNode->dfs_parent = &currentNode;
                this->dfs_stack.emplace_back(neighbor);
            }
        }
    }
    for (int i = (int) swaps.size() - 1; i >= 0; --i) {
        std::swap(neighborIds[i], neighborIds[swaps[i]]);
        swaps.pop_back();
    }
}

bool OuterplanarSubgraphDFS::CheckLeft(const NodeStruct &src,const NodeStruct &dst) {
    return !enclosingPointByLeftEdge(src, dst) && !crossesLeftDiagonal(dst);// || (dst.dfs_depth == src.last_spanned_root->dfs_depth - 1 && !current_left));
}

bool OuterplanarSubgraphDFS::CheckRight(const NodeStruct &src, const NodeStruct &dst) {
    return !enclosingPointByRightEdge(src, dst) && !crossesRightDiagonal(dst);// || (dst.dfs_depth == src.last_spanned_root->dfs_depth - 1 && !current_right));
}

bool OuterplanarSubgraphDFS::crossesLeftDiagonal(const NodeStruct& dst) {
    return !dst.nodeReachability.left();
}
bool OuterplanarSubgraphDFS::crossesRightDiagonal(const NodeStruct& dst) {
    return !dst.nodeReachability.right();
}

bool OuterplanarSubgraphDFS::enclosingPointByLeftEdge(const NodeStruct &src, const NodeStruct &dst) {
    return dst < src.last_L;
}

bool OuterplanarSubgraphDFS::enclosingPointByRightEdge(const NodeStruct &src, const NodeStruct &dst) {
    return dst < src.last_R;
}

void OuterplanarSubgraphDFS::GetValidEdges(NodeStruct& currentNode, std::vector<DirectedEdgeStruct>& edges, Side side) {
    edges.clear();
    for (NodeId neighborId : this->neighborsToVisit) {
        NodeStruct &neighbor = this->graphNodes[neighborId];
        if (side == Side::LEFT) {
            if (CheckLeft(currentNode, neighbor)) {
                edges.emplace_back(currentNode, neighbor);
            }
        } else {
            if (CheckRight(currentNode, neighbor)) {
                edges.emplace_back(currentNode, neighbor);
            }
        }
    }
}

void OuterplanarSubgraphDFS::AddEdges(const std::vector<DirectedEdgeStruct>& edges, Side side) {
    if (!edges.empty()) {
//        if (side == Side::LEFT){
//            left = true;
//        }
//        else{
//            right = true;
//        }

        for (auto &edge: edges) {
            NodeStruct& src = this->graphNodes[edge.src];
            NodeStruct& dst = this->graphNodes[edge.dst];
            //Update low
//            if (dst < *src.last_spanned_root){
//                left = true;
//                right = true;
//            }
            int i = (int) this->reachableNodes.size() - 2;
            NodeStruct *currentNode = &this->graphNodes[this->reachableNodes.back()];
            if (side == Side::LEFT) {
                currentNode->last_R = currentNode->dfs_parent->dfs_depth;
            } else {
                currentNode->last_L = currentNode->dfs_parent->dfs_depth;
            }
            while (!this->reachableNodes.empty() && this->graphNodes[this->reachableNodes[i]] > dst) {
                currentNode = &this->graphNodes[this->reachableNodes[i]];
                //Mark as side
                if (side == Side::LEFT) {
                    currentNode->nodeReachability.close_left();
                    currentNode->last_R = currentNode->dfs_depth;
                    currentNode->left_spanning = true;
                }
                else{
                    currentNode->nodeReachability.close_right();
                    currentNode->last_L = currentNode->dfs_depth;
                    currentNode->right_spanning = true;
                }
                this->reachableNodes.erase(this->reachableNodes.end() - 2);
                --i;
            }
        }
    }
}




void OuterplanarSubgraphDFS::reset() {
    for (int i = 0; i < this->_graph.nodes(); ++i) {
        this->graphNodes[i].reset();
    }
    dfs_stack.clear();
    reachableNodes.clear();
    neighborsToVisit.clear();
}

OuterplanarSubgraphDFS::OuterplanarSubgraphDFS(const OuterplanarSubgraphDFS& other) :  OuterplanarSubgraph(other) {
    reachableNodes = other.reachableNodes;
    graphNodes = other.graphNodes;
    dfs_root_node = other.dfs_root_node;
    possibleEdges = other.possibleEdges;
    neighborsToVisit = other.neighborsToVisit;
    dfs_stack = other.dfs_stack;
    left_edges = other.left_edges;
    right_edges = other.right_edges;
    //Used to take neighbor uniformly at random
    neighborIds = other.neighborIds;
    swaps = other.swaps;
    print = other.print;
    dfs_tree = other.dfs_tree;
}

void OuterplanarSubgraphDFS::UpdateNodeParameters(NodeStruct &currentNode) {
    if (print) {
        std::cout << " " << currentNode.dfs_parent->node_id << "-" << currentNode.node_id << " ";
    }
    //get parent node
    NodeStruct &currentParent = *currentNode.dfs_parent;

    //check if current node is in a new path
    if (this->graphNodes[this->reachableNodes.back()] > currentParent) {
        //Update the stack
        while (this->graphNodes[this->reachableNodes.back()] > currentParent) {
            this->reachableNodes.pop_back();
        }
        NodeStruct& path_root = currentParent;
        NodeStruct& path_root_parent = *path_root.dfs_parent;
        //if current branching node is not reachable from left (i.e. there is a left edge over the branching node) then there can be no higher left edge than to the branching get_node
        if (path_root.left_spanning) {
            path_root.last_L = path_root.dfs_depth;
        }
        else{
            path_root.last_L = path_root_parent.last_L;
        }
        //if current branching node is not reacheable from right (i.e there is a right edge over the branching node) then there can be no higher right edge than to the branching get_node
        if (path_root.right_spanning) {
            path_root.last_R = path_root.dfs_depth;
        }
        else{
            path_root.last_R = path_root_parent.last_R;
        }
        //Set reachability for branching get_node to true in branch
        path_root.nodeReachability.set(true, true);
        reachableNodes.emplace_back(path_root.node_id);
    }
    currentNode.last_L = currentParent.last_L;
    currentNode.last_R = currentParent.last_R;

    //Set get_node depth in dfs tree and id of the branch in which it was found
    currentNode.dfs_depth = currentParent.dfs_depth + 1;
}

void OuterplanarSubgraphDFS::AddGraphEdges(GraphStruct& subgraph, const std::vector<DirectedEdgeStruct> &edges, Side side) {
    for (auto &edge: edges) {
        NodeStruct src = this->graphNodes[edge.src];
        NodeStruct dst = this->graphNodes[edge.dst];
        if (print){
            if (side == Side::LEFT) {
                std::cout << " " << src.node_id << "--L" << dst.node_id << " ";
            }
            else{
                std::cout << " " << src.node_id << "--R" << dst.node_id << " ";
            }
        }
        subgraph.AddEdge(src.node_id, dst.node_id);
    }
}


#endif //CLOSURES_OUTERPLANARSUBGRAPHDFS_H