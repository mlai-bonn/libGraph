//
// Created by florian on 03.11.23.
//

#ifndef GRAPH_EXTENDED_H
#define GRAPH_EXTENDED_H

/**
 * @brief Extended GraphBase class
 * @details This class is used to extend the GraphBase class with additional functionality, it stores visited and distance vectors for faster graph search algorithms. User this class if you do multiple BFS or DFS searches on the same graph.
 * The drawback is that the memory consumption is higher as the visited and distance vectors and the parents are stored.
 */
class GraphExtended : public GraphStruct {
public:
    // default constructor
    GraphExtended() = default;
    GraphExtended(const std::string& name, INDEX size);
    GraphExtended(const std::string& name, NodeId size, const Labels& labels);
    explicit GraphExtended(const GraphStruct& graph);

    // read from path
    explicit GraphExtended(const std::string& path);

    void ResetSearchInformation();
public:
    NodeId _visited_id = 0;
    std::vector<NodeId> _visited;
    std::vector<NodeId> _distances;
    std::vector<NodeId> _parents;

};

inline void GraphExtended::ResetSearchInformation() {
    {
        ++this->_visited_id;
        if (this->_visited.size() != this->nodes()) {
            this->_visited_id = 0;
            this->_visited.resize(this->nodes(), this->_visited_id);
            std::fill(this->_visited.begin(), this->_visited.end(), this->_visited_id);
        }
        if (this->_visited_id == 0){
            std::fill(this->_visited.begin(), this->_visited.end(), this->_visited_id);
        }
        if (this->_distances.size() != this->nodes()) {
            this->_distances.resize(this->nodes(), -1);
        }
        std::fill(this->_distances.begin(), this->_distances.end(), -1);
    }
}

inline GraphExtended::GraphExtended(const std::string &path)  : GraphStruct(path) {

        _visited_id = 0;
        _visited = std::vector<NodeId>(this->nodes(), 0);
        _distances = std::vector<NodeId>(this->nodes(), 0);
        _parents = std::vector<NodeId>(this->nodes(), 0);
    }

inline GraphExtended::GraphExtended(const GraphStruct &graph)  : GraphStruct(graph) {
    _visited_id = 0;
    _visited = std::vector<NodeId>(graph.nodes(), 0);
    _distances = std::vector<NodeId>(graph.nodes(), 0);
    _parents = std::vector<NodeId>(graph.nodes(), 0);
}

inline GraphExtended::GraphExtended(const std::string& name, INDEX size) : GraphStruct(name, size) {
    _visited_id = 0;
    _visited = std::vector<NodeId>(size, 0);
    _distances = std::vector<NodeId>(size, 0);
    _parents = std::vector<NodeId>(size, 0);
}

inline GraphExtended::GraphExtended(const std::string& name, NodeId size, const Labels &labels) : GraphStruct(name, size, labels) {
    _visited_id = 0;
    _visited = std::vector<NodeId>(size, 0);
    _distances = std::vector<NodeId>(size, 0);
    _parents = std::vector<NodeId>(size, 0);
}
#endif //GRAPH_EXTENDED_H
