//
// Created by Florian on 15.04.2021.
//

#ifndef GRAPH_BASE_H
#define GRAPH_BASE_H
#include <vector>
#include <set>
#include <random>
#include <ranges>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "typedefs.h"
#include "GraphFunctions.h"
#include "LoadSave.h"


struct EdgeIterator;

/** Base class for graphs
 *
 */
struct GraphStruct{

    /**
     * Default constructor
     */
    GraphStruct()= default;

    /**
     * Default destructor
     */
    virtual ~GraphStruct() = default;

    /**
     * Constructor for creating a graph with the given name
     * @param name name of the graph
     * @param size number of nodes
     */
    explicit GraphStruct(const std::string& name, INDEX size);

    /**
     * Constructor for creating a graph with the given size and labels
     * @param name name of the graph
     * @param size number of nodes
     * @param labels labels of the nodes
     */
    GraphStruct(const std::string& name, NodeId size, const Labels& labels);

    /**
     * Constructor for loading a graph from a file
     * @param graphPath path to the graph file
     * @param relabeling if true the nodes will be relabeled from 0 to n-1
     * @param withLabels if true the graph will be loaded with labels
     * @param labelPath path to the label file
     * @param format format of the graph file
     * @param search_name if not empty only graphs with the given name will be loaded
     */
    explicit GraphStruct(const std::string & graphPath, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "", const std::string& format = "", const std::string& search_name = "");

    /**
     * Function for loading a graph from a file
     * @param graphPath path to the graph file
     * @param relabeling if true the nodes will be relabeled from 0 to n-1
     * @param withLabels if true the graph will be loaded with labels
     * @param labelPath path to the label file
     * @param format format of the graph file
     * @param search_name if not empty only graphs with the given name will be loaded
     */
    virtual void Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format = "", const std::string& search_name = "");

    /**
     * Function for saving the graph to a file
     * @param saveParams parameters for saving the graph
     */
    virtual void Save(const SaveParams& saveParams);

    /**
     * Function for initializing the graph
     * @param name
     * @param size
     * @param edges
     * @param nodeFeatures
     * @param edgeFeatures
     * @param nodeFeatureNames
     * @param edgeFeatureNames
     */
    virtual void Init(const std::string& name, int size, int edges, int nodeFeatures, int edgeFeatures, const std::vector<std::string>& nodeFeatureNames, const std::vector<std::string>& edgeFeatureNames);

    /**
     * Reset the graph to an empty graph with the given size
     * @param size
     */
    virtual void Reset(INDEX size);
    /**
     * Get the name of the graph
     * @return name of the graph
     */
    virtual std::string GetName() const{return _name;};

    // TODO check the following functions and comment them



    // Graph manipulation (functions that change the underlying graph data)
    // Nodes
    virtual INDEX AddNode();
    virtual INDEX AddNode(Label label);
    virtual INDEX AddNode(Label label, const std::vector<double>& nodeData);
    virtual INDEX AddNodes(INDEX number);
    virtual INDEX AddNodes(INDEX number, const std::vector<Label>& labels);
    virtual INDEX AddNodes(INDEX number, const std::vector<Label>& labels, const std::vector<std::vector<double>>& nodeData);

    virtual void RemoveNode(NodeId nodeId);
    virtual void RemoveNodes(const std::vector<NodeId>& nodeIds);

    virtual void RelabelNode(NodeId nodeId, Label newLabel);
    virtual void RelabelNode(NodeId nodeId, Label newLabel, const std::vector<double>& nodeData);

    // Edges
    virtual bool AddEdge(NodeId source, NodeId destination);
    virtual bool AddEdge(const std::pair<NodeId, NodeId>& edge);
    virtual bool AddEdge(NodeId source, NodeId destination, bool check_existence);
    virtual bool AddEdge(const std::pair<NodeId, NodeId>& edge, bool check_existence);
    virtual bool AddEdge(NodeId source, NodeId destination, const std::vector<double>& edgeData);
    virtual bool AddEdge(NodeId source, NodeId destination, const std::vector<double>& edgeData, bool check_existence);
    virtual bool AddEdge(const std::pair<NodeId, NodeId>& edge, const std::vector<double>& edgeData, bool check_existence);

    virtual bool RemoveEdge(NodeId source, NodeId destination);
    virtual bool RemoveEdge(const std::pair<NodeId, NodeId>& edge);

    virtual void RelabelEdge(NodeId source, NodeId destination, const std::vector<double>& edgeData){};
    virtual void RelabelEdge(const std::pair<NodeId, NodeId>& edge, std::vector<double>& newEdgeData){};

    /**
     * Resets the labels of the whole graph to the given labels
     * @param labels new labels for the graph
     */
    void SetLabels(const Labels *labels) ;

    /**
     * Set the name of the graph
     * @param name name of the graph
     */
    void SetName(const std::string& name) ;

    // Data access (const functions that do not change the underlying graph data)
    virtual double GetNodeData(NodeId node, const std::string& type) const {return 0;};
    virtual double GetNodeData(NodeId node, int index) const {return 0;};
    virtual const std::vector<double>& GetNodeData(NodeId node)const {return {0};};

    virtual double GetEdgeData(const EDGE & edge, const std::string& type) const {return 0;};
    virtual double GetEdgeData(const EDGE & edge, int index) const {return 0;};
    virtual const std::vector<double>& GetEdgeData(const EDGE& edge) const {return {0};};



    virtual void ReadNodeFeatures(double value, INDEX pos, const std::string& nodeFeatureName);

    virtual void WriteGraph(std::ofstream& Out, const SaveParams& saveParams);
    virtual void WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams);
    virtual void WriteEdges(std::ofstream& Out,const SaveParams& saveParams);

    bool GetConnectivity();




    // Setters and Getters
    INDEX GetNumLabels() const{return _labelSet.size();};

    void SetPath(const std::string& path){_path = path;};
    std::string GetPath(){return _path;};

    void SetType(GraphType type);

    GraphType GetType() const;

    virtual void write_graph_nodes(const std::string& graphPath, const std::string& fileName, const Nodes& nodes);


    const std::vector<Nodes>& graph() const {return _graph;};

    INDEX nodes() const;

    INDEX edges() const;

    void set_edge_num(INDEX edges){ _edges=edges;};

    virtual INDEX degree(NodeId nodeId) const;

    virtual INDEX degree(NodeId nodeId, Label label);

    virtual const std::vector<INDEX>& degreeByLabel(NodeId nodeId);

    virtual const std::vector<Label>& labels() const;

    virtual Label label(NodeId nodeId) const;

    virtual const Nodes& get_neighbors(NodeId nodeId) const;

    virtual Nodes& neighbors(NodeId nodeId);

    virtual NodeId neighbor(NodeId nodeId, INDEX neighborIdx) const;

    virtual NodeId neighbor(NodeId nodeId, INDEX neighborIdx, Label label);

    virtual NodeId random_neighbor_in_range(NodeId nodeId, INDEX minIdx, std::mt19937_64& gen);

    virtual void InitLabels();




    //Static functions
    //compare the given labeled degree vector with the labeled degree vector of the node with Id:srcId of this _graph
    bool isBigger(const std::vector<INDEX>& labeledDegree, NodeId nodeId);

    // overload the << operator for printing the graph
    friend std::ostream& operator<<(std::ostream& os, const GraphStruct& graph) {
        // print all the nodes with labels "Nodes: // 1, 0, ... "
        // print all edges between nodes, one edge per line  "Edges: // start_id -- edge_label[optional] -- [>] (optional if directed) end_id"
        // Print graph name
        os << "Graph: " << graph._name << std::endl;
        os << "\tType: " << graph.GetType() << std::endl;
        os << "\tLabel Type: " << graph.labelType << std::endl;
        os << "\t" << "Nodes: " << graph.nodes() << std::endl;
        // print nodes if there are labels
        if (graph.labelType != LABEL_TYPE::UNLABELED) {
            os << "\t" << "Labels: " << std::endl;
            // check wheter there are labels
            if (graph.labelType != LABEL_TYPE::UNLABELED) {
                os << "\t\t";
                for (INDEX i = 0; i < graph.nodes(); ++i) {
                    os << graph.label(i) << " ";
                }
                os << std::endl;
            }
        }
        os << "\t" << "Edges: " << graph.edges() << std::endl;
        // iterate over all nodes
        for (INDEX i = 0; i < graph.nodes(); ++i) {
            for (INDEX j = 0; j < graph.degree(i); ++j) {
                os << "\t\t";
                // print the node id
                os << i;
                // TODO print the edge label if given
                // print directed or undirected edge
                os << " -- ";
                // print the neighbor
                os << graph.neighbor(i, j) << std::endl;
            }
        }
        os << std::endl;
        return os;
    }


    //Iterators
    struct NodeIterator{
        // Prefix increment
        void operator++() { ++nodeId;};
        NodeId operator*() const { return nodeId; }
        friend bool operator== (const NodeIterator& a, INDEX b) { return a.nodeId == b;};
        friend bool operator!= (const NodeIterator& a, INDEX b) { return a.nodeId != b;};

        const GraphStruct* _graph{};
        NodeId nodeId{};
    };
    [[nodiscard]] NodeIterator begin() const { return NodeIterator{this, 0}; }
    [[nodiscard]] INDEX end() const {
        return this->nodes();
    }

    struct EdgeIterator{
        // Prefix increment
        void operator++() {
            ++neighborIdx;
            bool condition = (neighborIdx >= _graph->degree(srcId) || _graph->neighbor(srcId, neighborIdx) > srcId);
            while(condition){
                    ++srcId;
                    neighborIdx = 0;
                if (srcId >= this->_graph->nodes()){
                    return;
                }
                condition = (neighborIdx >= _graph->degree(srcId) || _graph->neighbor(srcId, neighborIdx) > srcId);
            }
            dstId = this->_graph->neighbor(srcId, neighborIdx);
        };

        std::pair<NodeId, NodeId> operator*() const {
            return {dstId, srcId};
        }
        friend bool operator== (const EdgeIterator& a, std::pair<NodeId, INDEX> b) {
            return a.srcId == b.first && a.neighborIdx == b.second;
        };
        friend bool operator!= (const EdgeIterator& a, std::pair<NodeId, INDEX> b) {
            return a.srcId != b.first || a.neighborIdx != b.second;
        };
        friend bool operator!= (const EdgeIterator& a, NodeId b) {
            return a.srcId != b;
        };


        const GraphStruct* _graph{};
        NodeId srcId{};
        NodeId dstId{};
        NodeId neighborIdx{};

    };

    [[nodiscard]] EdgeIterator first_edge() const {
        if (edges() == 0){
            return EdgeIterator(this, this->nodes(), 0,0);
        }
        NodeId startNode = 0;
        INDEX neighborIdx = 0;
        bool condition = neighborIdx >= this->degree(startNode) || this->neighbor(startNode, neighborIdx) > startNode;
        while(condition){
            ++startNode;
            neighborIdx = 0;
            if (startNode >= this->nodes()){
                break;
            }
            condition = neighborIdx >= this->degree(startNode) || this->neighbor(startNode, neighborIdx) > startNode;
        }
        return EdgeIterator{this,startNode, this->_graph[startNode][neighborIdx],  0};
    }

    [[nodiscard]] NodeId last_edge() const {
        return this->nodes();
    }

    struct EdgeIteratorDirected{
        // Prefix increment
        void operator++() {
            ++neighborIdx;
            bool condition = neighborIdx >= _graph->degree(srcId);
            while(condition){
                ++srcId;
                neighborIdx = 0;
                if (srcId >= this->_graph->nodes()){
                    return;
                }
                condition = neighborIdx >= _graph->degree(srcId);
            }
            dstId = this->_graph->neighbor(srcId, neighborIdx);
        };

        std::pair<NodeId, NodeId> operator*() const {
            return {dstId, srcId};
        }
        friend bool operator== (const EdgeIteratorDirected& a, std::pair<NodeId, INDEX> b) {
            return a.srcId == b.first && a.neighborIdx == b.second;
        };
        friend bool operator!= (const EdgeIteratorDirected& a, std::pair<NodeId, INDEX> b) {
            return a.srcId != b.first || a.neighborIdx != b.second;
        };
        friend bool operator!= (const EdgeIteratorDirected& a, NodeId b) {
            return a.srcId != b;
        };



        const GraphStruct* _graph{};
        NodeId srcId{};
        NodeId dstId{};
        NodeId neighborIdx{};

    };

    [[nodiscard]] EdgeIteratorDirected first_edge_directed() const {
        if (edges() == 0){
            return EdgeIteratorDirected(this, this->nodes(), 0,0);
        }
        NodeId startNode = 0;
        INDEX neighborIdx = 0;
        bool condition = neighborIdx >= this->degree(startNode);
        while(condition){
            ++startNode;
            neighborIdx = 0;
            if (startNode >= this->nodes()){
                break;
            }
            condition = neighborIdx >= this->degree(startNode);
        }
        return EdgeIteratorDirected{this,startNode, this->_graph[startNode][neighborIdx],  0};
    }

    [[nodiscard]] NodeId last_edge_directed() const {
        return this->nodes();
    }

    struct NeighborIterator{
        // Prefix increment
        void operator++() { ++_idx;};
        NodeId operator*() const {
            if (gen == nullptr) {
                return _graph->neighbor(_nodeId, _idx);
            }
            else{
                return _graph->random_neighbor_in_range(_nodeId, _idx, *gen);
            }
        }
        friend bool operator== (const NeighborIterator& a, INDEX b) { return a._idx == b;};
        friend bool operator!= (const NeighborIterator& a, INDEX b) { return a._idx != b;};

        GraphStruct* _graph;
        std::mt19937_64* gen;
        NodeId _nodeId;
        INDEX _idx;
    };

    NeighborIterator beginN(NodeId nodeId, std::mt19937_64* gen = nullptr) { return NeighborIterator{this, gen, nodeId, 0}; }
    [[nodiscard]] INDEX endN(NodeId nodeId, int maxIdx = -1) const   {
        if (maxIdx == -1) {
            return this->degree(nodeId);
        }
        else{
            return maxIdx;
        }
    }

    struct sort_by_name
    {
        inline bool operator() (const GraphStruct& struct1, const GraphStruct& struct2)
        {
            return (struct1._name < struct2._name);
        }
    };

    struct sort_by_size
    {
        inline bool operator() (const GraphStruct& struct1, const GraphStruct& struct2)
        {
            return (struct1.nodes() < struct2.nodes());
        }
    };

    bool operator==(const GraphStruct &rhs) const {
        if (_name != rhs._name){
            return false;
        }
        if (graphType != rhs.graphType){
            return false;
        }
        if (maxDegree != rhs.maxDegree){
            return false;
        }
        if (labelType != rhs.labelType){
            return false;
        }
        if (GetNumLabels() != rhs.GetNumLabels()){
            return false;
        }
        if (labelMap != rhs.labelMap){
            return false;
        }
        if (labelFrequencyMap != rhs.labelFrequencyMap){
            return false;
        }
        if (_nodes != rhs._nodes){
            return false;
        }
        if (_edges != rhs._edges){
            return false;
        }
        if (_graph != rhs._graph){
            return false;
        }
        if (_degrees != rhs._degrees){
            return false;
        }
        if (_labels != rhs._labels){
            return false;
        }
        if (labeledDegreeMap != rhs.labeledDegreeMap){
            return false;
        }
        if (labeledDegreeVector != rhs.labeledDegreeVector){
            return false;
        }
        if (labeledVectorGraph != rhs.labeledVectorGraph){
            return false;
        }
        if (labeledMapGraph != rhs.labeledMapGraph){
            return false;
        }
        return true;
    }

public:
    //Graph attributes for faster access
    int maxDegree = 0;
    LABEL_TYPE labelType = LABEL_TYPE::UNLABELED;

    //for labels
    std::unordered_map<Label, Nodes> labelMap{};
    std::unordered_map<Label, INDEX> labelFrequencyMap{};
    std::pair<std::unordered_map<Label, Label>, std::unordered_map<Label, Label>> labelMappingPair{};

    //store original node Ids
    std::unordered_map<INDEX ,INDEX> IdsToOriginalIds;

    virtual bool IsEdge(NodeId source, NodeId destination) const;
    virtual bool IsEdge(const std::pair<NodeId, NodeId>& edge) const;
    virtual bool IsEdge(const EdgeIterator& edge) const;
    virtual bool AddEdge(const EdgeIterator& edge);
    virtual bool AddEdge(const EdgeIterator& edge, bool check_existence);
    virtual bool RemoveEdge(const EdgeIterator& edge);

    //Get get_neighbors by []
    const Nodes& operator[](const NodeId nodeId){return _graph[nodeId];};
    void sortNeighborIds();
    bool has_neighbor_label(NodeId nodeId, Label label);

    bool comp_degree(NodeId i, NodeId j) const{
        return this->degree(i) < this->degree(j);
    };

    void degree_sort(Nodes& nodes) const {
        std::ranges::sort(nodes, [this](const NodeId l, const NodeId r){return comp_degree(l,r);});
    }







protected:
    // General variables
    std::string _name;
    std::string _path;
    GraphType graphType = GraphType::GENERAL;

    // Changing under graph manipulation
    INDEX _nodes = 0;
    INDEX _edges = 0;
    std::vector<Nodes> _graph;
    std::vector<INDEX> _degrees;

    bool _checked_connected = false;
    bool _is_connected = false;



    //for labels
    Labels _labels = Labels();
    std::set<Label> _labelSet;
    std::vector<std::unordered_map<Label, INDEX>> labeledDegreeMap{};
    std::vector<std::vector<INDEX>> labeledDegreeVector{};
    std::vector<std::vector<Nodes>> labeledVectorGraph{};
    std::vector<std::unordered_map<Label, Nodes>> labeledMapGraph{};

    // protected functions


    /**
     * Creates a new graph from name size and labels (TODO merge with Init)
     * @param name name of the graph
     * @param size number of nodes
     * @param labels labels of the nodes
     */
    void create_graph(const std::string& name, NodeId size, const Labels& labels);

    /**
     * Updates the graph structure after adding or removing nodes or edges
     */
    virtual void update_graph_struct();

    /**
     * Updates the graph node labels after adding or removing nodes or relabeling nodes
     */
    void update_node_label_information();



    void set_graph(const std::vector<Nodes>& nodes);


public :
    [[deprecated]]
    void Save_0(const std::string & graphPath = "", GraphFormat Format = GraphFormat::BINARY, bool Labeled = false, bool OnlyGraph = true, const std::string& Name = "") const;



};






inline GraphStruct::GraphStruct(const std::string& name, const INDEX size) {
    create_graph(name, size, Labels());
}

/// Construct new undirected _graph
/// \param size
/// \param labels
inline GraphStruct::GraphStruct(const std::string& name, const INDEX size, const Labels& labels) {
    create_graph(name, size, labels);
}

/// Initialize a new undirected _graph
/// \param size
/// \param labels
inline void GraphStruct::create_graph(const std::string& name, INDEX size, const Labels& labels) {
    this->SetName(name);
    this->_nodes = size;
    this->_edges = 0;
    this->SetLabels(&labels);
    this->_graph = std::vector<std::vector<INDEX>>(nodes(), std::vector<INDEX>());
    this->_degrees = std::vector<INDEX>(nodes(), 0);
}


inline void GraphStruct::update_graph_struct() {
    this->update_node_label_information();
}

/// Initialize labeled _graph using the label type
/// \param label
inline void GraphStruct::update_node_label_information() {
    // save original labels in mapping
    this->labelMap = GraphFunctions::GetGraphLabelMap(_labels);
    this->labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
    this->labelMappingPair = GraphFunctions::GetLabelMappingPair(_labels);
    this->_labelSet = std::set<Label>{};
    for (auto const& label : _labels) {
        this->_labelSet.insert(label);
    }
    // Map original label to new label
    //this->labelMappingPair.first = std::unordered_map<Label, Label>();
    //for (Label label = 0; label < this->_labels.size(); ++label) {}
    // map new label to original label
    //this->labelMappingPair.second = std::unordered_map<Label, Label>();
    //for (auto const& [label, index] : this->labelMap) {
    //    this->labelMappingPair.second[index] = label;
    //}


    // If there are more than 10 labels use dense representation otherwise sparse (map vs vector)
    if (GetNumLabels() >= 10){
        labelType = LABEL_TYPE::LABELED_DENSE;
        this->labeledMapGraph.resize(nodes(),std::unordered_map<Label, Nodes>{});
        this->labeledDegreeMap.resize(nodes(), std::unordered_map<Label, INDEX>{});
    }
    else{
        labelType = LABEL_TYPE::LABELED_SPARSE;
        this->labeledVectorGraph.resize( nodes(), std::vector<std::vector<INDEX>>());
        for (NodeId Id = 0; Id < nodes(); ++Id) {
            this->labeledVectorGraph[Id].resize(GetNumLabels(), std::vector<INDEX>{});
        }
        this->labeledDegreeVector.resize(nodes(), std::vector<INDEX>());
        for (NodeId Id = 0; Id < nodes(); ++Id) {
            this->labeledDegreeVector[Id].resize(GetNumLabels(), 0);
        }
    }

    //set labels of neighbor Nodes
    for (NodeId Id = 0; Id < nodes(); ++Id) {
        for (NodeId neighborId: _graph[Id]) {
            //Create labeled _graph
            if (labelType == LABEL_TYPE::LABELED_DENSE) {
                Label NodeLabel = _labels[neighborId];
                if (this->labeledMapGraph[Id].find(NodeLabel)== this->labeledMapGraph[Id].end()) {
                    Nodes newNodes = Nodes();
                    newNodes.push_back(neighborId);
                    this->labeledMapGraph[Id][NodeLabel] = newNodes;
                    this->labeledDegreeMap[Id][NodeLabel] = 1;
                } else {
                    this->labeledMapGraph[Id][NodeLabel].emplace_back(neighborId);
                    ++this->labeledDegreeMap[Id][NodeLabel];
                }
            } else if (labelType == LABEL_TYPE::LABELED_SPARSE) {
                Label NodeLabel = _labels[neighborId];
                // Map node label to id using the label map
                Label MappedNodeLabel = this->labelMappingPair.first[NodeLabel];
                this->labeledVectorGraph[Id][MappedNodeLabel].emplace_back(neighborId);
                ++this->labeledDegreeVector[Id][MappedNodeLabel];
            }
        }
    }
}

/// Get the degree of a _graph node
/// \param nodeId
/// \return
inline INDEX GraphStruct::degree(NodeId nodeId) const {
    return this->get_neighbors(nodeId).size();
    //return this->_degrees[srcId];
}

/// Get the number of _graph nodes
/// \return
inline INDEX GraphStruct::nodes() const {
    return _nodes;
}

/// Get the number of _graph edges
/// \return
inline INDEX GraphStruct::edges() const {
    return _edges;
}



/// Check if a _graph contains some edge
/// \param source
/// \param destination
/// \return
inline bool GraphStruct::IsEdge(NodeId source, NodeId destination) const {
    if (this->degree(source) < this->degree(destination)){
        return std::binary_search(_graph[source].begin(), _graph[source].end(), destination);
    }
    else{
        return std::binary_search(_graph[destination].begin(), _graph[destination].end(), source);
    }
}

/**
 * Check if a _graph contains some _edge
 * \param edge
 */
inline bool GraphStruct::IsEdge(const std::pair<NodeId, NodeId>& edge) const {
    return IsEdge(edge.first, edge.second);
}

/**
 * Check if a _graph contains some _edge
 * \param edge
 */
inline bool GraphStruct::IsEdge(const EdgeIterator& edge) const {
    return IsEdge(*edge);
}

inline bool GraphStruct::AddEdge(const EdgeIterator &edge) {
    return AddEdge(*edge, true);
}

/// Get all neighbors of a _graph node (const)
/// \param nodeId
/// \return
inline const Nodes &GraphStruct::get_neighbors(const NodeId nodeId) const {
    return _graph[nodeId];
}

/// Get all neighbors of a _graph node (non const)
/// \param nodeId
/// \return
inline Nodes &GraphStruct::neighbors(const NodeId nodeId) {
    return _graph[nodeId];
}

/// Add an edge to an undirected _graph TODO give computational complexity
/// \param source
/// \param destination
/// \param check_existence
/// \return
inline bool GraphStruct::AddEdge(NodeId source, NodeId destination, const bool check_existence) {
    if (!check_existence || !IsEdge(source, destination)){
        this->_graph[source].emplace_back(destination);
        INDEX ElementId = static_cast<INDEX>(this->_graph[source].size()) - 1;
        while (ElementId > 0 && this->_graph[source][ElementId] < this->_graph[source][ElementId - 1]){
            std::swap(this->_graph[source][ElementId], this->_graph[source][ElementId - 1]);
            --ElementId;
        }

        this->_graph[destination].emplace_back(source);
        ElementId = static_cast<INDEX>(this->_graph[destination].size()) - 1;
        while (ElementId > 0 && this->_graph[destination][ElementId] < this->_graph[destination][ElementId - 1]){
            std::swap(this->_graph[destination][ElementId], this->_graph[destination][ElementId - 1]);
            --ElementId;
        }
        this->_checked_connected = false;
        ++this->_edges;
        ++this->_degrees[source];
        ++this->_degrees[destination];
        this->maxDegree = std::max(this->maxDegree, static_cast<int>(std::max(this->degree(source), this->degree(destination))));
        return true;
    }
    return false;
}

/// Add an edge to an undirected _graph
/// \param edge
/// \param check_existence
/// \param source
/// \param destination
/// \return
inline bool GraphStruct::AddEdge(const std::pair<NodeId, NodeId>& edge, const bool check_existence) {
    return AddEdge(edge.first, edge.second, check_existence);
}

inline bool GraphStruct::AddEdge(const NodeId source, const NodeId destination, const std::vector<double> &edgeData) {
    return AddEdge(source, destination, true);
}

inline bool GraphStruct::AddEdge(const NodeId source, const NodeId destination, const std::vector<double>& edgeData,
                                 const bool check_existence) {
    return AddEdge(source, destination, check_existence);
}

inline bool GraphStruct::AddEdge(const std::pair<NodeId, NodeId> &edge, const std::vector<double>& edgeData,
    const bool check_existence) {
    return AddEdge(edge.first, edge.second, edgeData, check_existence);
}

/// Add an edge to an undirected _graph
/// \param edge
/// \param check_existence
/// \return
inline bool GraphStruct::AddEdge(const EdgeIterator& edge, const bool check_existence) {
    return AddEdge(*edge, check_existence);
}

inline void GraphStruct::set_graph(const std::vector<Nodes> &nodes) {
    this->_nodes = nodes.size();
    this->_edges = 0;
    this->_graph = nodes;
    this->_degrees.resize(_nodes);
    for (NodeId i = 0; i < _nodes; ++i) {
        const Nodes& neighbors = nodes[i];
        this->_degrees[i] = neighbors.size();
        this->_edges += neighbors.size();
    }
    this->_edges /= 2;
    this->sortNeighborIds();
}

/// Get labels of the _graph nodes
/// \return
inline const Labels& GraphStruct::labels() const {
    return _labels;
}

/// Get label of a _graph node
/// \param nodeId
/// \return
inline Label GraphStruct::label(NodeId nodeId) const{
    return _labels[nodeId];
}


/// Get specific neighbor of a _graph node
/// \param nodeId
/// \param neighborIdx
/// \return
inline NodeId GraphStruct::neighbor(NodeId nodeId, INDEX neighborIdx) const {
    return _graph[nodeId][neighborIdx];
}

/// Get specific neighbor of a _graph node with specific label
/// \param nodeId
/// \param neighborIdx
/// \param label
/// \return
inline NodeId GraphStruct::neighbor(NodeId nodeId, INDEX neighborIdx, Label label) {
    if (labelType == LABEL_TYPE::UNLABELED){
        return neighbor(nodeId, neighborIdx);
    }
    else if (labelType == LABEL_TYPE::LABELED_SPARSE){
        return this->labeledVectorGraph[nodeId][label][neighborIdx];
    }
    else if(labelType == LABEL_TYPE::LABELED_DENSE){
        return this->labeledMapGraph[nodeId][label][neighborIdx];
    }
    return -1;
}

/// Compare the labeled degrees of _graph node with given labeled degrees
/// \param labeledDegree
/// \param nodeId
/// \return
inline bool GraphStruct::isBigger(const std::vector<INDEX>& labeledDegree, NodeId nodeId) {
    for (Label label = 0; label < labeledDegree.size(); ++label) {
        if (labeledDegreeVector[nodeId][label] > labeledDegree[label]){
            return true;
        }
    }
    return false;
}

/// sort all the neighbors according to their Ids
inline void GraphStruct::sortNeighborIds() {
    for (Nodes& nodes1 : this->_graph) {
        std::sort(nodes1.begin(), nodes1.end());
    }
    for (auto& [label, nodes1] : this->labelMap){
        std::sort(nodes1.begin(), nodes1.end());
    }
    for (std::unordered_map<Label, Nodes>& nodeLabelMap : this->labeledMapGraph){
        for (auto& [label, nodes1] : nodeLabelMap) {
            std::sort(nodes1.begin(), nodes1.end());
        }
    }
    for(std::vector<Nodes>& nodeLabelVector : this->labeledVectorGraph){
        for(Nodes& nodes1 : nodeLabelVector){
            std::sort(nodes1.begin(), nodes1.end());
        }
    }
}

/// Get labeled degree of specific node (all neighbors with specific label)
/// \param nodeId
/// \param label
/// \return
inline INDEX GraphStruct::degree(NodeId nodeId, Label label) {
    if (labelType == LABEL_TYPE::LABELED_DENSE){
        return this->labeledDegreeMap[nodeId][label];
    }
    else if (labelType == LABEL_TYPE::LABELED_SPARSE){
        return this->labeledDegreeVector[nodeId][label];
    }
    return -1;
}

/// Get degrees of a node specified by labels of the neighbors
/// \param nodeId
/// \return
inline const std::vector<INDEX> &GraphStruct::degreeByLabel(NodeId nodeId) {
    return this->labeledDegreeVector[nodeId];
}


/// Check if _graph node has neighbor with specific label
/// \param nodeId
/// \param label
/// \return
inline bool GraphStruct::has_neighbor_label(NodeId nodeId, Label label) {
    if (labelType == LABEL_TYPE::LABELED_SPARSE){
        return this->labeledVectorGraph[nodeId].size() > label;
    }
    else if(labelType == LABEL_TYPE::LABELED_DENSE){
        return this->labeledMapGraph[nodeId].find(label) != this->labeledMapGraph[nodeId].end();
    }
    return false;
}

/// Get random neighbor of a node in a specific range (changes the order of the nodes)
/// \param nodeId
/// \param minIdx
/// \param gen
/// \return
inline NodeId GraphStruct::random_neighbor_in_range(NodeId nodeId, INDEX minIdx, std::mt19937_64& gen) {
    INDEX randIdx = std::uniform_int_distribution<INDEX>(minIdx, this->degree(nodeId) - 1)(gen);
    std::swap(this->neighbors(nodeId)[randIdx], this->neighbors(nodeId)[0]);
    return this->neighbors(nodeId)[0];
}

/// Construct _graph from path (optionally with node labels)
/// \param graphPath
/// \param relabeling
/// \param withLabels
/// \param labelPath
inline GraphStruct::GraphStruct(const std::string &graphPath, bool relabeling, bool withLabels, const std::string & labelPath, const std::string& format, const std::string& search_name) {
    GraphStruct::Load(graphPath, relabeling, withLabels, labelPath, format, search_name);
}

/// Save _graph in certain path and format
/// \param graphPath
/// \param Format
/// \param Labeled
/// \param OnlyGraph
/// \param Name
inline void GraphStruct::Save(const SaveParams& saveParams) {
    std::string saveName;
    if (saveParams.graphPath.empty()) {
        if (!this->_path.empty()) {
            saveName = this->_path;
        } else {
            return;
        }
    }
    else {
        saveName = saveParams.graphPath;
    }
    if (!saveParams.Name.empty()) {
        saveName += saveParams.Name;
    }
    else{
        saveName += this->_name;
    }
    switch (saveParams.Format) {
        case GraphFormat::BGF: {
            int SaveVersion = 1;
            int graphNumber = 1;
            std::ofstream Out(saveName + ".bgf", std::ios::out | std::ios::binary);
            Out.write((char *) (&SaveVersion), sizeof(int));
            Out.write((char *) (&graphNumber), sizeof(int));
            for (int i = 0; i < graphNumber; ++i) {
                this->WriteGraph(Out, saveParams);
            }

            for (int i = 0; i < graphNumber; ++i) {
                this->WriteNodeFeatures(Out, saveParams);
                this->WriteEdges(Out, saveParams);
            }
            Out.close();
            break;
        }
        case GraphFormat::BGFS: {
            int SaveVersion = 1;
            int graphNumber = 1;
            std::ofstream Out(saveName + ".bgfs", std::ios::out | std::ios::binary);
            Out.write((char *) (&SaveVersion), sizeof(int));
            Out.write((char *) (&graphNumber), sizeof(int));
            for (int i = 0; i < graphNumber; ++i) {
                this->WriteGraph(Out, saveParams);
            }

            for (int i = 0; i < graphNumber; ++i) {
                this->WriteNodeFeatures(Out, saveParams);
                this->WriteEdges(Out, saveParams);
            }
            Out.close();
            break;
        }
        case GraphFormat::BINARY: {
            std::ofstream Out(saveName + ".bin", std::ios::out | std::ios::binary);
            const std::string graphName = this->_name;
            unsigned int stringLength = graphName.length();
            GraphType Type = this->graphType;
            Out.write((char *) (&stringLength), sizeof(stringLength));
            Out.write(graphName.c_str(), stringLength);
            Out.write((char *) (&Type), sizeof(GraphType));
            auto Size = (INDEX) this->nodes();
            INDEX Edges = this->edges();
            Out.write((char *) (&Size), sizeof(INDEX));
            Out.write((char *) (&Edges), sizeof(INDEX));
            INDEX Src = 0;
            for (auto const &edges: this->_graph) {
                for (auto Dst: edges) {
                    if (Dst > Src) {
                        Out.write((char *) (&Src), sizeof(INDEX));
                        Out.write((char *) (&Dst), sizeof(INDEX));
                    }
                }
                ++Src;
            }
            Out.write((char *) (&saveParams.Labeled), sizeof(bool));
            if (saveParams.Labeled) {
                for (auto const &L: labels()) {
                    Out.write((char *) (&L), sizeof(Label));
                }
            }
            Out.write((char *) (true), sizeof(bool));
            Out.close();
            break;
        }
        case GraphFormat::EDGES: {
            std::ofstream Out(saveName + ".edges", std::ios::out);
            INDEX Src = 0;
            bool FirstLine = true;
            for (auto const &edges: this->_graph) {
                for (auto Dst: edges) {
                    if (Dst > Src) {
                        if (!FirstLine) {
                            Out << "\n";
                        }
                        Out << Src << " " << Dst;
                        FirstLine = false;
                    }
                }
                ++Src;
            }
            break;
        }
        case GraphFormat::PEREGRINE_SMALL: {
            std::cout << saveName + ".peregrine";
            std::ofstream Out(saveName + ".peregrine", std::ios::out);
            INDEX id = 0;
            if (labelType != LABEL_TYPE::UNLABELED) {
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    for (auto Node: Nodes) {
                        if (NodeCounter < Node) {
                            if (EdgeCounter > 0) {
                                Out << std::endl;
                            }
                            Out << NodeCounter << " " << label(NodeCounter) << " " << Node << " " << label(Node);
                            ++EdgeCounter;
                        }
                    }
                    ++NodeCounter;
                }
            } else {
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    for (auto Node: Nodes) {
                        if (EdgeCounter > 0) {
                            Out << std::endl;
                        }
                        if (NodeCounter < Node) {
                            Out << NodeCounter << " " << Node;
                        }
                        ++EdgeCounter;
                    }
                    ++NodeCounter;
                }

            }
            Out.close();
            break;
        }
        case GraphFormat::PEREGRINE_DATA: {
            std::cout << saveName + ".peregrine" << std::endl;
            std::ofstream Out(saveName + ".peregrine", std::ios::out);
            INDEX id = 0;
            if (labelType != LABEL_TYPE::UNLABELED) {
                std::ofstream OutLabels(saveParams.graphPath + this->_name + ".peregrine_label", std::ios::out);
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    OutLabels << NodeCounter << label(NodeCounter) << std::endl;
                    for (auto Node: Nodes) {
                        if (NodeCounter < Node) {
                            if (EdgeCounter > 0) {
                                Out << std::endl;
                            }
                            Out << NodeCounter << " " << Node;
                            ++EdgeCounter;
                        }
                    }
                    ++NodeCounter;
                }
            } else {
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    for (auto Node: Nodes) {
                        if (NodeCounter < Node) {
                            if (EdgeCounter > 0) {
                                Out << std::endl;
                            }
                            Out << NodeCounter << " " << Node;
                            //std::cout << NodeCounter << " " << Node << std::endl;
                            ++EdgeCounter;
                        }
                    }
                    ++NodeCounter;
                }
            }
            Out.close();
            break;
        }
        case GraphFormat::DIMACS:
            break;
        case GraphFormat::AIDS:
            break;
    }
}



inline void GraphStruct::ReadNodeFeatures(double value, INDEX pos, const std::string &nodeFeatureName) {
    if (nodeFeatureName == "label") {
        this->_labels.emplace_back(value);
    }
}

inline void GraphStruct::Init(const std::string &name, int size, int edges, int nodeFeatures, int edgeFeatures,
                              const std::vector<std::string> &nodeFeatureNames,
                              const std::vector<std::string> &edgeFeatureNames) {
    this->_name = name;
    this->_nodes = size;
    this->_edges = edges;
    this->_degrees.resize(_nodes);
    this->_graph.resize(_nodes);

}

/// Write a vector of nodes to a file
/// \param graphPath
/// \param fileName
/// \param nodes
inline void GraphStruct::write_graph_nodes(const std::string & graphPath, const std::string& fileName, const Nodes& nodes) {
    std::ofstream Out(graphPath + fileName + ".nodes", std::ios::out);
    bool FirstLine = true;
    for (NodeId x : nodes) {
        if (FirstLine) {
            Out << this->IdsToOriginalIds[x];
            FirstLine = false;
        }
        else{
            Out << " " << this->IdsToOriginalIds[x];
        }
    }
    Out.close();
}

/// Add labels to a _graph
/// \param labels
inline void GraphStruct::SetLabels(const Labels* labels) {
    // return an error if the number of nodes does not coincide with the labels
    if (labels->size() == this->_graph.size()){
        this->_labels = *labels;
        update_node_label_information();
    }
    else{
        //TODO throw exception
        std::cerr << "Error: The number of labels does not coincide with the number of nodes." << std::endl;
    }
}

inline void GraphStruct::Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format, const std::string& search_name) {
    int graphId = 0;
    std::string path = graphPath;
    if (!search_name.empty()) {
        std::vector<std::string> possible_files;
        std::vector<std::string> files;
        for (const auto &entry: std::filesystem::directory_iterator(path)) {
            if (entry.is_regular_file()) {
                possible_files.emplace_back(entry.path().string());
            }
        }
        for (const auto &p: possible_files) {
            std::string lower_name = search_name;
            std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                           [](unsigned char c) { return std::tolower(c); } // correct
            );
            std::string path_lower = p;
            std::transform(path_lower.begin(), path_lower.end(), path_lower.begin(),
                           [](unsigned char c) { return std::tolower(c); } // correct
            );
            if (path_lower.find(lower_name) != std::string::npos) {
                files.emplace_back(p);
            }
        }
        if (files.size() == 1){
            path = files[0];
        }
    }

    if (std::filesystem::is_regular_file(path)) {
        this->_path = std::filesystem::path(path).remove_filename().string();
        std::string extension = std::filesystem::path(path).extension().string();
        if (extension == ".bgf" || extension == ".bgfs") {
            int saveVersion = 1;
            int graphNumber;
            std::vector<std::string> graphsNames;
            std::vector<GraphType> graphsTypes;
            std::vector<INDEX> graphsSizes;
            std::vector<std::vector<std::string>> graphsNodeFeatureNames;
            std::vector<INDEX> graphsEdges;
            std::vector<std::vector<std::string>> graphsEdgeFeatureNames;
            std::ifstream In(path, std::ios::in | std::ios::binary);

            if (LoadSave::ReadBGF(extension, In, saveVersion, graphNumber, graphsNames, graphsTypes, graphsSizes,
                                     graphsNodeFeatureNames, graphsEdges, graphsEdgeFeatureNames)) {
                for (int i = 0; i < graphNumber; ++i) {
                    //Create _graph
                    if (i == graphId) {
                        this->Init(graphsNames[i], (int) graphsSizes[i], (int) graphsEdges[i],
                                   (int) graphsNodeFeatureNames[i].size(), (int) graphsEdgeFeatureNames[i].size(),
                                   graphsNodeFeatureNames[i], graphsEdgeFeatureNames[i]);
                    }
                    //Read the nodes
                    for (int j = 0; j < graphsSizes[i]; ++j) {
                        for (int k = 0; k < graphsNodeFeatureNames[i].size(); ++k) {
                            double val;
                            if (extension == ".bgf") {
                                In.read((char *) (&val), sizeof(double));
                            } else if (extension == ".bgfs") {
                                unsigned int int_val;
                                In.read((char *) (&int_val), sizeof(unsigned int));
                                val = int_val;
                            }
                            if (i == graphId) {
                                this->ReadNodeFeatures(val, j, graphsNodeFeatureNames[i][k]);
                            }
                        }
                    }
                    //Read the edges
                    INDEX added_edges = 0;
                    INDEX original_edges = this->_edges;
                    for (INDEX j = 0; j < original_edges; ++j) {
                        INDEX Src = 0;
                        INDEX Dst = 0;
                        if (extension == ".bgf") {
                            In.read((char *) (&Src), sizeof(INDEX));
                            In.read((char *) (&Dst), sizeof(INDEX));
                        } else if (extension == ".bgfs") {
                            unsigned int int_Src;
                            unsigned int int_Dst;
                            In.read((char *) (&int_Src), sizeof(unsigned int));
                            In.read((char *) (&int_Dst), sizeof(unsigned int));
                            Src = int_Src;
                            Dst = int_Dst;
                        }

                        std::vector<double> edgeData;
                        for (int k = 0; k < graphsEdgeFeatureNames[i].size(); ++k) {
                            double val;
                            In.read((char *) (&val), sizeof(double));
                            edgeData.emplace_back(val);
                        }
                        if (graphId == i) {
                            if (this->AddEdge(Src, Dst, edgeData)) {
                                ++added_edges;
                            }
                        }
                    }
                    this->_edges = added_edges;

                    if (graphId == i) {
                        this->InitLabels();
                    }
                }
            }
            In.close();
        }
        else if (extension == ".bin") {
            std::ifstream In(path, std::ios::in | std::ios::binary);
            std::string Name;
            unsigned int stringLength;
            In.read((char *) (&stringLength), sizeof(stringLength));
            Name.resize(stringLength);
            In.read((char *) Name.c_str(), stringLength);
            this->_name = Name;
            In.read((char *) (&this->graphType), sizeof(GraphType));
            In.read((char *) (&this->_nodes), sizeof(INDEX));
            In.read((char *) (&this->_edges), sizeof(INDEX));
            for (INDEX i = 0; i < nodes(); ++i) {
                this->_graph.emplace_back();
            }
            for (INDEX i = 0; i < edges(); ++i) {
                INDEX Src = 0;
                INDEX Dst = 0;
                In.read((char *) (&Src), sizeof(INDEX));
                In.read((char *) (&Dst), sizeof(INDEX));
                if (this->_graph.size() > std::max(Src, Dst)) {
                    this->_graph[Src].emplace_back(Dst);
                    this->_graph[Dst].emplace_back(Src);
                } else {
                    INDEX x = std::max(Src, Dst);
                }
            }
            this->_degrees = std::vector<INDEX>(nodes(), 0);
            INDEX Id = 0;
            for (auto const &graphNode: this->_graph) {
                this->_degrees[Id] = (INDEX) graphNode.size();
                if (this->_degrees[Id] > this->maxDegree) {
                    this->maxDegree = static_cast<int>(this->_degrees[Id]);
                }
                ++Id;
            }
            bool Labeled;
            In.read((char *) (&Labeled), sizeof(bool));
            if (Labeled && withLabels) {
                Label L;
                for (INDEX i = 0; i < nodes(); ++i) {
                    In.read((char *) (&L), sizeof(Label));
                    this->_labels.emplace_back(L);
                }
                if (this->_labels.size() == this->_graph.size()) {
                    update_node_label_information();
                } else {
                    //TODO throw exception
                }

            }
            bool OnlyGraph;
            In.read((char *) (&OnlyGraph), sizeof(bool));
            if (!OnlyGraph) {
            }
            In.close();
        } else if (extension == ".edges" || extension == ".txt") {
            std::string graph_name = std::filesystem::path(path).stem().string();
            this->_name = graph_name;
            NodeId src;
            NodeId dest;
            std::string a, b;
            std::string line;
            std::ifstream infile(path);
            std::vector<std::pair<INDEX, INDEX>> graphEdges;
            std::set<INDEX> graphNodeIds;
            std::unordered_map<INDEX, INDEX> originalIdsToNodeIds;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                iss >> a;
                if (a == "#") {
                    continue;
                }
                iss >> b;
                // Check if file is of format
                // num_nodes
                // node_idA node_idB
                // ...
                if (b.empty()) {
                    for (INDEX i = 0; i < std::stoull(a); ++i) {
                        this->_graph.emplace_back();
                    }
                } else {
                    src = std::stoull(a);
                    dest = std::stoull(b);
                    graphEdges.emplace_back(src, dest);
                    graphNodeIds.emplace(src);
                    graphNodeIds.emplace(dest);
                }
            }
            _nodes = graphNodeIds.size();
            NodeId num_edges = graphEdges.size();
            this->_degrees.resize(_nodes);
            this->_graph.resize(_nodes);
            INDEX nodeCounter = 0;
            for (auto x: graphNodeIds) {
                originalIdsToNodeIds.insert({x, nodeCounter});
                ++nodeCounter;
            }
            unsigned int num_edges_duplicates = 0;
            for (auto edge: graphEdges) {
                if (!GraphStruct::AddEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second], true)) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (this->IsEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base _graph)!"
                      << std::endl;
            std::cout << graphEdges.size() - num_edges_duplicates << " edges are remaining" << std::endl;
            if (!labelPath.empty() && withLabels) {
                std::string label_extension = std::filesystem::path(labelPath).extension().string();
                std::ifstream label_file(labelPath);
                if (label_extension == ".vertexids") {
                    std::string label_line;
                    _labels = Labels(nodes(), 0);
                    Label label = 0;
                    while (std::getline(label_file, label_line)) {
                        std::istringstream iss(label_line);
                        std::string token;
                        while (std::getline(iss, token, ' ')) {
                            _labels[originalIdsToNodeIds[std::stoi(token)]] = label;
                        }
                        ++label;
                    }
                } else {
                    try {
                        NodeId id;
                        Label label;
                        while (label_file >> id >> label) {
                            _labels.push_back(label);
                        }
                    }
                    catch (...) {
                    }
                }
                if (this->_labels.size() == this->_graph.size()) {
                    update_node_label_information();
                } else {
                    //TODO throw exception
                }
            }
        } else if (format == "dimacs") {
            std::string graph_name = std::filesystem::path(path).stem().string();
            this->_name = graph_name;
            NodeId src;
            NodeId dest;
            std::string a, b;
            std::string line;
            std::ifstream infile(path);
            std::vector<std::pair<INDEX, INDEX>> graphEdges;
            std::set<INDEX> graphNodeIds;
            std::unordered_map<INDEX, INDEX> originalIdsToNodeIds;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                iss >> a;
                if (a == "p") {
                    continue;
                }
                iss >> a;
                iss >> b;
                // Check if file is of format
                // num_nodes
                // node_idA node_idB
                // ...
                if (b.empty()) {
                    for (INDEX i = 0; i < std::stoull(a); ++i) {
                        this->_graph.emplace_back();
                    }
                } else {
                    src = std::stoull(a);
                    dest = std::stoull(b);
                    graphEdges.emplace_back(src, dest);
                    graphNodeIds.emplace(src);
                    graphNodeIds.emplace(dest);
                }
            }
            _nodes = graphNodeIds.size();
            NodeId num_edges = graphEdges.size();
            this->_degrees.resize(_nodes);
            this->_graph.resize(_nodes);
            INDEX nodeCounter = 0;
            for (auto x: graphNodeIds) {
                originalIdsToNodeIds.insert({x, nodeCounter});
                ++nodeCounter;
            }
            int num_edges_duplicates = 0;
            for (auto edge: graphEdges) {
                if (!GraphStruct::AddEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second], true)) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (this->IsEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            if (num_edges_duplicates > 0) {
                std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base _graph)!"
                          << std::endl;
            }
            if (!labelPath.empty() && withLabels) {
                std::string label_extension = std::filesystem::path(labelPath).extension().string();
                std::ifstream label_file(labelPath);
                if (label_extension == ".vertexids") {
                    std::string label_line;
                    _labels = Labels(nodes(), 0);
                    Label label = 0;
                    while (std::getline(label_file, label_line)) {
                        std::istringstream iss(label_line);
                        std::string token;
                        while (std::getline(iss, token, ' ')) {
                            _labels[originalIdsToNodeIds[std::stoi(token)]] = label;
                        }
                        ++label;
                    }
                } else {
                    try {
                        NodeId id;
                        Label label;
                        while (label_file >> id >> label) {
                            _labels.push_back(label);
                        }
                    }
                    catch (...) {
                    }
                }
                if (this->_labels.size() == this->_graph.size()) {
                    update_node_label_information();
                } else {
                    //TODO throw exception
                }
            }
        }else if (format == "aids") {
            std::string graph_name = std::filesystem::path(path).stem().string();
            this->_name = graph_name;
            NodeId src;
            NodeId dest;
            std::string a, b;
            std::string line;
            std::ifstream infile(path);
            std::vector<std::pair<INDEX, INDEX>> graphEdges;
            std::set<INDEX> graphNodeIds;
            std::unordered_map<INDEX, INDEX> originalIdsToNodeIds;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                iss >> a;
                if (a == "p") {
                    continue;
                }
                iss >> a;
                iss >> b;
                // Check if file is of format
                // num_nodes
                // node_idA node_idB
                // ...
                if (b.empty()) {
                    for (INDEX i = 0; i < std::stoull(a); ++i) {
                        this->_graph.emplace_back();
                    }
                } else {
                    src = std::stoull(a);
                    dest = std::stoull(b);
                    graphEdges.emplace_back(src, dest);
                    graphNodeIds.emplace(src);
                    graphNodeIds.emplace(dest);
                }
            }
            _nodes = graphNodeIds.size();
            NodeId num_edges = graphEdges.size();
            this->_degrees.resize(_nodes);
            this->_graph.resize(_nodes);
            INDEX nodeCounter = 0;
            for (auto x: graphNodeIds) {
                originalIdsToNodeIds.insert({x, nodeCounter});
                ++nodeCounter;
            }
            int num_edges_duplicates = 0;
            for (auto edge: graphEdges) {
                if (!GraphStruct::AddEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second], true)) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (this->IsEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            if (num_edges_duplicates > 0) {
                std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base _graph)!"
                          << std::endl;
            }
            if (!labelPath.empty() && withLabels) {
                std::string label_extension = std::filesystem::path(labelPath).extension().string();
                std::ifstream label_file(labelPath);
                if (label_extension == ".vertexids") {
                    std::string label_line;
                    _labels = Labels(nodes(), 0);
                    Label label = 0;
                    while (std::getline(label_file, label_line)) {
                        std::istringstream iss(label_line);
                        std::string token;
                        while (std::getline(iss, token, ' ')) {
                            _labels[originalIdsToNodeIds[std::stoi(token)]] = label;
                        }
                        ++label;
                    }
                } else {
                    try {
                        NodeId id;
                        Label label;
                        while (label_file >> id >> label) {
                            _labels.push_back(label);
                        }
                    }
                    catch (...) {
                    }
                }
                if (this->_labels.size() == this->_graph.size()) {
                    update_node_label_information();
                } else {
                    //TODO throw exception
                }
            }
        }
    }

    this->sortNeighborIds();
    if (GraphFunctions::CheckTree(*this)){
        this->SetType(GraphType::TREE);
    }

}

/// Deprecated Save _graph in certain path and format
/// \param graphPath
/// \param Format
/// \param Labeled
/// \param OnlyGraph
/// \param Name
inline void GraphStruct::Save_0(const std::string &graphPath, GraphFormat Format, bool Labeled, bool OnlyGraph, const std::string& Name) const {
    int SaveVersion = 0;
    std::string saveName;
    if (graphPath.empty()) {
        if (!this->_path.empty()) {
            saveName = this->_path + this->_name;
        } else {
            return;
        }
    }
    else {
        saveName = graphPath + this->_name;
    }
    if (!Name.empty()) {
        saveName = graphPath + Name;
    }
    switch (Format) {
        case GraphFormat::BINARY: {
            std::ofstream Out(saveName + ".bin", std::ios::out | std::ios::binary);
            const std::string graphName = this->_name;
            unsigned int stringLength = graphName.length();
            GraphType Type = this->graphType;
            Out.write((char *) (&stringLength), sizeof(stringLength));
            Out.write(graphName.c_str(), stringLength);
            Out.write((char *) (&Type), sizeof(GraphType));
            auto Size = (INDEX) this->nodes();
            INDEX Edges = this->edges();
            Out.write((char *) (&Size), sizeof(INDEX));
            Out.write((char *) (&Edges), sizeof(INDEX));
            INDEX Src = 0;
            for (auto const &edges: this->_graph) {
                for (auto Dst: edges) {
                    if (Dst > Src) {
                        Out.write((char *) (&Src), sizeof(INDEX));
                        Out.write((char *) (&Dst), sizeof(INDEX));
                    }
                }
                ++Src;
            }
            Out.write((char *) (&Labeled), sizeof(bool));
            if (Labeled) {
                for (auto const &L: labels()) {
                    Out.write((char *) (&L), sizeof(Label));
                }
            }
            Out.write((char *) (&OnlyGraph), sizeof(bool));
            if (!OnlyGraph) {
            }
            Out.close();
            break;
        }
        case GraphFormat::EDGES: {
            std::ofstream Out(saveName + ".edges", std::ios::out);
            INDEX Src = 0;
            bool FirstLine = true;
            for (auto const &edges: this->_graph) {
                for (auto Dst: edges) {
                    if (Dst > Src) {
                        if (!FirstLine) {
                            Out << "\n";
                        }
                        Out << Src << " " << Dst;
                        FirstLine = false;
                    }
                }
                ++Src;
            }
        }
        case GraphFormat::PEREGRINE_SMALL: {
            std::cout << saveName + ".peregrine";
            std::ofstream Out(saveName + ".peregrine", std::ios::out);
            INDEX id = 0;
            if (labelType != LABEL_TYPE::UNLABELED) {
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    for (auto Node: Nodes) {
                        if (NodeCounter < Node) {
                            if (EdgeCounter > 0) {
                                Out << std::endl;
                            }
                            Out << NodeCounter << " " << label(NodeCounter) << " " << Node << " " << label(Node);
                            ++EdgeCounter;
                        }
                    }
                    ++NodeCounter;
                }
            } else {
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    for (auto Node: Nodes) {
                        if (EdgeCounter > 0) {
                            Out << std::endl;
                        }
                        if (NodeCounter < Node) {
                            Out << NodeCounter << " " << Node;
                        }
                        ++EdgeCounter;
                    }
                    ++NodeCounter;
                }

            }
            Out.close();
            break;
        }
        case GraphFormat::PEREGRINE_DATA: {
            std::cout << saveName + ".peregrine" << std::endl;
            std::ofstream Out(saveName + ".peregrine", std::ios::out);
            INDEX id = 0;
            if (labelType != LABEL_TYPE::UNLABELED) {
                std::ofstream OutLabels(graphPath + this->_name + ".peregrine_label", std::ios::out);
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    OutLabels << NodeCounter << label(NodeCounter) << std::endl;
                    for (auto Node: Nodes) {
                        if (NodeCounter < Node) {
                            if (EdgeCounter > 0) {
                                Out << std::endl;
                            }
                            Out << NodeCounter << " " << Node;
                            ++EdgeCounter;
                        }
                    }
                    ++NodeCounter;
                }
            } else {
                NodeId NodeCounter = 0;
                INDEX EdgeCounter = 0;
                for (auto const &Nodes: this->_graph) {
                    for (auto Node: Nodes) {
                        if (NodeCounter < Node) {
                            if (EdgeCounter > 0) {
                                Out << std::endl;
                            }
                            Out << NodeCounter << " " << Node;
                            //std::cout << NodeCounter << " " << Node << std::endl;
                            ++EdgeCounter;
                        }
                    }
                    ++NodeCounter;
                }
            }
            Out.close();
            break;
        }
    }
}

inline void GraphStruct::SetName(const std::string &name) {
    this->_name = name;
}



inline void GraphStruct::InitLabels() {
    //Add _graph labels if possible
    if (this->_labels.size() == this->_graph.size()) {
        update_node_label_information();
    } else {
        this->labelType = LABEL_TYPE::UNLABELED;
    }
}

inline void GraphStruct::WriteGraph(std::ofstream& Out, const SaveParams& saveParams) {
    const std::string graphName = this->_name;
    unsigned int stringLength = graphName.length();
    GraphType Type = this->graphType;
    auto Size = this->nodes();
    unsigned int nodeFeatures;
    std::vector<std::string> nodeFeatureNames;

    Out.write(reinterpret_cast<char *>(&stringLength), sizeof(unsigned int));
    Out.write(graphName.c_str(), stringLength);
    Out.write(reinterpret_cast<char *>(&Type), sizeof(GraphType));

    if(saveParams.Format == GraphFormat::BGF) {
        Out.write(reinterpret_cast<char *>(&Size), sizeof(INDEX));
    }
    else if(saveParams.Format == GraphFormat::BGFS) {
        unsigned int int_size = static_cast<unsigned int>(Size);
        Out.write(reinterpret_cast<char *>(&int_size), sizeof(int_size));
    }

    nodeFeatures = 0;
    if (saveParams.Labeled || _labels.size() == _nodes) {
        nodeFeatures = 1;
        nodeFeatureNames.emplace_back("label");
    }
    Out.write(reinterpret_cast<char *>(&nodeFeatures), sizeof(unsigned int));

    for (int j = 0; j < nodeFeatures; ++j) {
        unsigned int nodeFeatureStringLength = nodeFeatureNames[j].length();
        Out.write(reinterpret_cast<char *>(&nodeFeatureStringLength), sizeof(unsigned int));
        Out.write(nodeFeatureNames[j].c_str(), nodeFeatureStringLength);
    }

    auto Edges = this->edges();
    if(saveParams.Format == GraphFormat::BGF) {
        Out.write(reinterpret_cast<char *>(&Edges), sizeof(INDEX));
    }
    else if(saveParams.Format == GraphFormat::BGFS) {
        auto int_edges = (unsigned int) Edges;
        Out.write(reinterpret_cast<char *>(&int_edges), sizeof(unsigned int));
    }

    unsigned int edgeFeatures = 0;
    std::vector<std::string> edgeFeatureNames = {};
    Out.write(reinterpret_cast<char *>(&edgeFeatures), sizeof(unsigned int));
    for (int j = 0; j < edgeFeatures; ++j) {
        unsigned int edgeFeatureStringLength = edgeFeatureNames[j].length();
        Out.write(reinterpret_cast<char *>(&edgeFeatureStringLength), sizeof(unsigned int));
        Out.write(edgeFeatureNames[j].c_str(), edgeFeatureStringLength);
    }

}

inline void GraphStruct::WriteNodeFeatures(std::ofstream& Out, const SaveParams &saveParams) {
    unsigned int nodeFeatures = 0;
    if (saveParams.Labeled || this->_labels.size() == nodes()) {
        nodeFeatures = 1;
        for (int j = 0; j < this->nodes(); ++j) {
            for (int k = 0; k < nodeFeatures; ++k) {
                if (saveParams.Format == GraphFormat::BGF) {
                    auto label = (double) this->labels()[j];
                    Out.write((char *) (&label), sizeof(double));
                }
                else if(saveParams.Format == GraphFormat::BGFS){
                    auto label = (unsigned int) this->labels()[j];
                    Out.write((char *) (&label), sizeof(unsigned int));
                }
            }
        }
    }

}

inline void GraphStruct::WriteEdges(std::ofstream& Out, const SaveParams &saveParams) {
    INDEX Src = 0;
    for (auto const &edges: this->_graph) {
        for (auto Dst: edges) {
            if (Dst > Src) {
                if (saveParams.Format == GraphFormat::BGF) {
                    Out.write((char *) (&Src), sizeof(INDEX));
                    Out.write((char *) (&Dst), sizeof(INDEX));
                }
                else if (saveParams.Format == GraphFormat::BGFS){
                    auto int_src = static_cast<unsigned int>(Src);
                    auto int_dst = static_cast<unsigned int>(Dst);
                    Out.write((char *) (&int_src), sizeof(unsigned int));
                    Out.write((char *) (&int_dst), sizeof(unsigned int));
                }
                unsigned int edgeFeatures = 0;
                for (int j = 0; j < edgeFeatures; ++j) {
                    //Out.write((char *) (&Src), sizeof(INDEX));
                }
            }
        }
        ++Src;
    }
}

inline bool GraphStruct::GetConnectivity() {
    if (this->_checked_connected) {
        return this->_is_connected;
    }

    std::vector<INDEX> reachability = std::vector<INDEX>(this->nodes(), -1);
    INDEX number;
    this->_is_connected = GraphFunctions::ReachableNodes(*this, 0, reachability, 0, number);
    this->_checked_connected = true;

    return this->_is_connected;
}

inline void GraphStruct::SetType(GraphType type) {
    this->graphType = type;
}

inline  GraphType GraphStruct::GetType() const {
    return this->graphType;
}

inline INDEX GraphStruct::AddNodes(const INDEX number) {
    // Emplace back the graph adjacency list
    for(INDEX i = 0; i < number; i++) {
        this->_graph.emplace_back();
        this->_degrees.emplace_back(0);
        ++this->_nodes;
        this->_labels.emplace_back(0); // default label 0
    }
    // if nodes are added the _graph is not connected anymore
    if (number > 0) {
        this->_is_connected = false;
    }
    this->update_graph_struct();
    return static_cast<INDEX>(this->_graph.size()) - 1;
}

inline INDEX GraphStruct::AddNodes(const INDEX number, const std::vector<Label>& labels) {
    GraphStruct::AddNodes(number);
    // If the label vector matches the number of nodes
    if (labels.size() == number) {
        for (INDEX i = 0; i < number; ++i) {
            this->_labels[this->_nodes - number + i] = labels[i];
        }
    }
    // node data is only used in the child class GraphLabeledBase
    this->update_graph_struct();
    return static_cast<INDEX>(this->_graph.size()) - 1;
}

inline INDEX GraphStruct::AddNodes(const INDEX number, const std::vector<Label>& labels,
    const std::vector<std::vector<double>>& nodeData) {
    return GraphStruct::AddNodes(number, labels);
}

inline INDEX GraphStruct::AddNode() {
    return GraphStruct::AddNodes(1);
}

inline INDEX GraphStruct::AddNode(const Label label) {
    return GraphStruct::AddNodes(1, std::vector<Label>{label});
}

inline INDEX GraphStruct::AddNode(const Label label, const std::vector<double>& nodeData) {
    return GraphStruct::AddNodes(1, std::vector<Label>{label});
}

inline void GraphStruct::RemoveNode(const NodeId nodeId) {
    // Raise error if node does not have degree 0
    if (this->_degrees[nodeId] != 0) {
        std::cerr << "Error: Node " << nodeId << " has degree " << this->_degrees[nodeId] << " and cannot be removed!" << std::endl;
        return;
    }
    --this->_nodes;
    // Iterate over all nodes and remove nodeId and decrease Id by 1 if Id is greater than nodeId
    INDEX node_counter = 0;
    for (auto& nodes : this->_graph) {
        bool found = false;
        for (INDEX j = 0; j < nodes.size(); ++j) {
            if (nodes[j] == nodeId) {
                found = true;
            }
            else {
                if (nodes[j] > nodeId) {
                    --nodes[j];
                }
                if (found) {
                    nodes[j-1] = nodes[j];
                }
            }
        }
        if (found) {
            nodes.pop_back();
            --this->_degrees[node_counter];
        }
        ++node_counter;
    }
    this->_graph.erase(this->_graph.begin() + nodeId);
    this->_degrees.erase(this->_degrees.begin() + nodeId);
    this->_labels.erase(this->_labels.begin() + nodeId);
    this->update_graph_struct();
}

inline void GraphStruct::RemoveNodes(const std::vector<NodeId> &nodeIds) {
    for (const auto nodeId : nodeIds) {
        this->RemoveNode(nodeId);
    }
}

inline void GraphStruct::RelabelNode(const NodeId nodeId, const Label newLabel) {
    this->_labels[nodeId] = newLabel;
    // TODO make the update more efficient (based on update event)
    this->update_node_label_information();
}

inline void GraphStruct::RelabelNode(const NodeId nodeId, const Label newLabel, const std::vector<double>& nodeData) {
    this->RelabelNode(nodeId, newLabel);
}

inline bool GraphStruct::AddEdge(const NodeId source, const NodeId destination) {
    return GraphStruct::AddEdge(source, destination, true);
}

inline bool GraphStruct::AddEdge(const std::pair<NodeId, NodeId> &edge) {
    return GraphStruct::AddEdge(edge.first, edge.second, true);
}


inline bool GraphStruct::RemoveEdge(NodeId source, NodeId destination) {
    auto it_source = std::find(this->_graph[source].begin(), this->_graph[source].end(), destination);
    auto it_destination = std::find(this->_graph[destination].begin(), this->_graph[destination].end(), source);
    if (it_source != this->_graph[source].end() && it_destination != this->_graph[destination].end()) {
        this->_graph[source].erase(it_source);
        this->_graph[destination].erase(it_destination);
        --_degrees[source];
        --_degrees[destination];
        --this->_edges;
        this->_checked_connected = false;
        return true;
    }
    return false;
}

inline bool GraphStruct::RemoveEdge(const std::pair<NodeId, NodeId> &edge) {
    return GraphStruct::RemoveEdge(edge.first, edge.second);
}

inline bool GraphStruct::RemoveEdge(const EdgeIterator &edge) {
    return GraphStruct::RemoveEdge(*edge);
}

inline void GraphStruct::Reset(const INDEX size) {
    this->_graph.resize(size);
    this->_degrees.resize(size);
    this->_nodes = size;
    this->_edges = 0;

    for (int i = 0; i < this->nodes(); ++i) {
        this->_graph[i].clear();
        this->_degrees[i] = 0;
    }
}







#endif //GRAPH_BASE_H
