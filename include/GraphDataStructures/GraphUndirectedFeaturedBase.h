//
// Created by florian on 23.10.23.
//

#ifndef TESTGRAPHLIB_GRAPHLABELEDBASE_H
#define TESTGRAPHLIB_GRAPHLABELEDBASE_H
#include "GraphUndirectedBase.h"

struct UDataGraph : UGraphStruct{
    /**
     * Default constructor
     */
    UDataGraph();
    /**
     * Default destructor
     */
    ~UDataGraph() override = default;

    /**
     *
     * @param name name of the graph
     * @param size size of the graph, i.e., number of nodes
     */
    UDataGraph(const std::string& name, INDEX size);

    /**
     *
     * @param name name of the graph
     * @param size size of the graph, i.e., number of nodes
     * @param labels labels of the nodes
     */
    UDataGraph(const std::string& name, INDEX size, const Labels& labels);

    /**
     *
     * @param graphPath
     * @param relabeling
     * @param withLabels
     * @param labelPath
     */
    explicit UDataGraph(const std::string & graphPath, const std::string& labelPath = "", bool relabeling = true, bool withLabels = false);


    // TODO check functions
    void Load(const std::string & graphPath, bool relabeling, bool withLabels, const std::string& labelPath, const std::string& format = "");
    void Save(const SaveParams& saveParams) override;

    void Init(const std::string& name, int size, int edges, int nodeFeatures, int edgeFeatures, const std::vector<std::string>& nodeFeatureNames, const std::vector<std::string>& edgeFeatureNames) override;
    void ReadNodeFeatures(double value, INDEX pos, const std::string& nodeFeatureName) override;

    void WriteGraph(std::ofstream& Out, const SaveParams& saveParams) override;
    void WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams) override;
    void WriteEdges(std::ofstream& Out,const SaveParams& saveParams) override;

    // Graph Manipulation
    // Graph manipulation (functions that change the underlying graph data)
    // Nodes
    INDEX AddNodes(INDEX number) override;

    INDEX AddNodes(INDEX number, const std::vector<Label>& labels) override;
    INDEX AddNodes(INDEX number, const std::vector<Label>& labels, const std::vector<std::vector<double>>& nodeData) override;

    void RemoveNode(NodeId nodeId) override;
    void RelabelNode(NodeId nodeId, Label newLabel) override;
    void RelabelNode(NodeId nodeId, Label newLabel, const std::vector<double>& nodeData) override;

    // Edges
    bool AddEdge(NodeId source, NodeId destination, const std::vector<double>& edgeData) override;
    bool AddEdge(NodeId source, NodeId destination, const std::vector<double>& edgeData, bool check_existence) override;
    bool AddEdge(const std::pair<NodeId, NodeId>& edge, const std::vector<double>& edgeData, bool check_existence) override;

    bool RemoveEdge(NodeId source, NodeId destination) override;
    bool RemoveEdge(const std::pair<NodeId, NodeId>& edge) override;

    void RelabelEdge(NodeId source, NodeId destination, const std::vector<double>& edgeData) override;
    void RelabelEdge(const std::pair<NodeId, NodeId>& edge, std::vector<double>& newEdgeData) override;

    // Data access (const functions that do not change the underlying graph data)
    double GetNodeData(NodeId node, const std::string& type) const override;
    double GetNodeData(NodeId node, int index) const override;
    const std::vector<double>& GetNodeData(NodeId node) const override;

    double GetEdgeData(const EDGE & edge, const std::string& type) const override;
    double GetEdgeData(const EDGE & edge, int index) const override;
    const std::vector<double>& GetEdgeData(const EDGE& edge) const override;


    void add_edge_data(const EDGE & edge, const std::string& type, double data);
    void add_edge_data(const EDGE & edge, int index, double data);
    void add_edge_data(const EDGE & edge, const std::vector<double>& data);

    // set node data names
    void set_node_data_names(const std::unordered_map<std::string, int>& nodeDataNames);
    // set edge data names
    void set_edge_data_names(const std::unordered_map<std::string, int>& edgeDataNames);


    static bool dijkstra(const UDataGraph& graph, NodeId src, double (*weight_function)(const UDataGraph&, EDGE), std::vector<double>& distances);
    static bool dijkstra(const UDataGraph& graph, NodeId src, NodeId dest, double (*weight_function)(const UDataGraph&, EDGE), std::vector<NodeId>& path, std::vector<double>& distances, double& length);

    bool operator==(const UDataGraph &rhs) const {
        if ((UGraphStruct) *this != (UGraphStruct) rhs){
            return false;
        }
        if (_edge_data_size != rhs._edge_data_size){
            return false;
        }
        if (_node_data_size != rhs._node_data_size){
            return false;
        }
        if (_node_data_names != rhs._node_data_names){
            return false;
        }
        if (_edge_data_names != rhs._edge_data_names){
            return false;
        }
        if (_node_data != rhs._node_data){
            return false;
        }
        if (_edge_data != rhs._edge_data){
            return false;
        }
        return true;
    }

private:
    std::unordered_map<std::string, int> _node_data_names;
    std::unordered_map<std::string, int> _edge_data_names;
    std::vector<std::vector<double>> _node_data;
    std::unordered_map<INDEX, std::unordered_map<INDEX, std::vector<double>>> _edge_data;
    int _edge_data_size = 0;
    int _node_data_size = 0;
};




/// Default constructor
inline UDataGraph::UDataGraph() = default;

inline UDataGraph::UDataGraph(const std::string &name, const INDEX size) : UGraphStruct(name, size) {
}

inline UDataGraph::UDataGraph(const std::string& name, const INDEX size, const Labels &labels) : UGraphStruct(name, size, labels) {
}

/// Read directed data _graph from file TODO binary extension
/// \param graphPath path of the _graph file
/// \param relabeling if true node labels start with 0 and end with size-1, otherwise use original labels TODO check this
/// \param withLabels if considering node attributes
/// \param labelPath path to node data
inline UDataGraph::UDataGraph(const std::string &graphPath, const std::string &labelPath, bool relabeling, bool withLabels){
    UDataGraph::Load(graphPath, relabeling, withLabels,labelPath);
}

inline void UDataGraph::Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format) {
    int Version = 1;
    int graphId = 0;
    if (std::filesystem::is_regular_file(graphPath)) {
        std::string extension = std::filesystem::path(graphPath).extension().string();
        if (extension == ".bgf"){
            UGraphStruct::Load(graphPath, relabeling,withLabels,labelPath);
        }
        else if (extension == ".edges" || extension == ".txt") {
            std::string graph_name = std::filesystem::path(graphPath).stem().string();
            this->_name = graph_name;
            NodeId src;
            NodeId dest;
            std::string a;
            std::string line;
            std::ifstream infile(graphPath);
            std::vector<std::pair<INDEX, INDEX>> graphEdges;
            std::set<INDEX> graphNodeIds;
            std::unordered_map<INDEX, std::unordered_map<INDEX, std::vector<double>>> original_edge_data;
            std::unordered_map<INDEX, INDEX> originalIdsToNodeIds;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                int count = 0;
                for (std::string x; iss >> x;) {
                    if (x == "#") {
                        int ct = 0;
                        for (; iss >> a;) {
                            if(ct > 1) {
                                this->_edge_data_names.insert({a, this->_edge_data_size});
                                ++this->_edge_data_size;
                            }
                            ++ct;
                        }
                        break;
                    } else {
                        if (count == 0) {
                            src = std::stoull(x);
                        } else if (count == 1) {
                            dest = std::stoull(x);
                            graphEdges.emplace_back(src, dest);
                            graphNodeIds.emplace(src);
                            graphNodeIds.emplace(dest);
                        } else {
                            original_edge_data[src][dest].emplace_back(std::stoull(x));
                        }
                        ++count;
                    }
                }
            }
            _nodes = graphNodeIds.size();
            this->_degrees.resize(_nodes);
            this->_graph.resize(_nodes);
            INDEX nodeCounter = 0;
            for (auto x: graphNodeIds) {
                originalIdsToNodeIds.insert({x, nodeCounter});
                this->IdsToOriginalIds.insert({nodeCounter, x});
                ++nodeCounter;
            }
            for (auto edge: graphEdges) {
                AddEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second], original_edge_data[edge.first][edge.second]);
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
                } else if (label_extension == ".data") {
                    try {
                        std::ifstream labelFile(labelPath);
                        NodeId id;
                        Label label;
                        this->_node_data.resize(_nodes);
                        while (std::getline(labelFile, line)) {
                            std::istringstream iss(line);
                            int count = 0;
                            for (std::string x; iss >> x;) {
                                if (x == "#") {
                                    int ct = 0;
                                    for (; iss >> a;) {
                                        if (ct > 0){
                                            this->_node_data_names.insert({a, this->_node_data_size});
                                            ++this->_node_data_size;
                                        }
                                        ++ct;
                                        ++count;
                                    }
                                    break;
                                } else {
                                    if (count == 0) {
                                        id = std::stoull(x);
                                    } else {
                                        if (this->_node_data_names.find("label") != this->_node_data_names.end() &&
                                            this->_node_data_names["label"] == count - 1) {
                                            this->_labels.resize(_nodes);
                                            this->_labels[originalIdsToNodeIds[id]] = std::stoull(x);
                                        }
                                        this->_node_data[originalIdsToNodeIds[id]].emplace_back(std::stoull(x));
                                    }
                                    ++count;
                                }
                            }
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

        this->sortNeighborIds();
        if (CheckTree()){
            this->graphType = GraphType::TREE;
        }
    }
    else{
        std::cout << "Graph path is not a valid path!" << std::endl;
    }
}

inline void UDataGraph::Save(const SaveParams& saveParams) {
    UGraphStruct::Save(saveParams);
}

/// Get the data of an edge using the data type
/// \param edge
/// \param type
/// \return
inline double UDataGraph::GetEdgeData(const EDGE & edge, const std::string& type) const {
    const int index = _edge_data_names.at(type);
    if (index == -1) {
        return 0;
    }
    return _edge_data.at(edge.first).at(edge.second)[index];
}

/// Get the data of an edge using the data index
/// \param edge
/// \param index
/// \return
inline double UDataGraph::GetEdgeData(const EDGE &edge, const int index) const {
    return _edge_data.at(edge.first).at(edge.second)[index];
}

/// Get edge data for some edge in a directed data _graph
/// \param edge
/// \return
inline const std::vector<double>& UDataGraph::GetEdgeData(const EDGE & edge) const {
    return _edge_data.at(edge.first).at(edge.second);
}


/// Add the data of an edge using the type
/// \param edge
/// \param type
/// \param data
inline void UDataGraph::add_edge_data(const EDGE &edge, const std::string &type, double data) {
    add_edge_data(edge, _edge_data_names[type], data);
}

/// Add the data of an edge using the edge index
/// \param edge
/// \param index
/// \param data
inline void UDataGraph::add_edge_data(const EDGE &edge, int index, double data) {
    _edge_data[edge.first][edge.second].resize(this->_edge_data_size);
    _edge_data[edge.first][edge.second][index] = data;
    _edge_data[edge.second][edge.first].resize(this->_edge_data_size);
    _edge_data[edge.second][edge.first][index] = data;
}

/// Add the data of an edge
/// \param edge
/// \param data
inline void UDataGraph::add_edge_data(const EDGE &edge, const std::vector<double> &data) {
    _edge_data[edge.first][edge.second] = data;
    _edge_data[edge.second][edge.first] = data;
}


inline void UDataGraph::set_node_data_names(const std::unordered_map<std::string, int> &nodeDataNames) {
    this->_node_data_names = nodeDataNames;
    this->_node_data_size = (int) nodeDataNames.size();
}

inline void UDataGraph::set_edge_data_names(const std::unordered_map<std::string, int> &edgeDataNames) {
    this->_edge_data_names = edgeDataNames;
    this->_edge_data_size = (int) edgeDataNames.size();
}

inline INDEX UDataGraph::AddNodes(const INDEX number, const std::vector<Label>& labels,
    const std::vector<std::vector<double>>& nodeData) {
    const INDEX result =  UGraphStruct::AddNodes(number, labels, nodeData);
    // Check whether nodeData size equals the number and add the data
    if ( nodeData.size() == number) {
        for (INDEX i = 0; i < number; ++i) {
            this->_node_data.emplace_back(nodeData.at(i));
        }
    }
    return result;
}

inline void UDataGraph::RelabelNode(const NodeId nodeId, const Label newLabel, const std::vector<double>& nodeData) {
    this->_node_data[nodeId] = nodeData;
    GraphStruct::RelabelNode(nodeId, newLabel, nodeData);
}

inline bool UDataGraph::AddEdge(const NodeId source, const NodeId destination, const std::vector<double> &edgeData) {
    return UDataGraph::AddEdge(source, destination, edgeData, true);
}

inline bool UDataGraph::AddEdge(NodeId source, NodeId destination, const std::vector<double>& edgeData,
                                bool check_existence) {
    add_edge_data(EDGE{source, destination}, edgeData);
    return UGraphStruct::AddEdge(source, destination);

}

inline bool UDataGraph::AddEdge(const std::pair<NodeId, NodeId> &edge, const std::vector<double> &edgeData,
    const bool check_existence) {
    return UDataGraph::AddEdge(edge.first, edge.second, edgeData, check_existence);
}

inline bool UDataGraph::RemoveEdge(const NodeId source, const NodeId destination) {
    this->_edge_data[source].erase(destination);
    this->_edge_data[destination].erase(source);
    return UGraphStruct::RemoveEdge(source, destination);
}

inline bool UDataGraph::RemoveEdge(const std::pair<NodeId, NodeId> &edge) {
    return UDataGraph::RemoveEdge(edge.first, edge.second);
}

inline void UDataGraph::RelabelEdge(const NodeId source, const NodeId destination, const std::vector<double> &edgeData) {
    this->_edge_data[source][destination] = edgeData;
    UGraphStruct::RelabelEdge(source, destination, edgeData);
}

inline void UDataGraph::RelabelEdge(const std::pair<NodeId, NodeId> &edge, std::vector<double> &newEdgeData) {
    UDataGraph::RelabelEdge(edge.first, edge.second, newEdgeData);
}

inline void UDataGraph::RemoveNode(const NodeId nodeId) {
    UGraphStruct::RemoveNode(nodeId);
    this->_node_data.erase(this->_node_data.begin() + nodeId);
    this->_edge_data.erase(nodeId);
    std::unordered_map<INDEX, std::unordered_map<INDEX, std::vector<double>>> new_edge_data;
    // Update the edge data iterate over all edge data
    for (auto & edge_map : this->_edge_data) {
        NodeId source = edge_map.first;
        std::unordered_map<INDEX, std::vector<double>> new_target_data_map;
        for (auto & [target, data] : edge_map.second) {
            if (target >= nodeId) {
                new_target_data_map[target - 1] = data;
            }
            else{
                new_target_data_map[target] = data;
            }
        }
        if (source >= nodeId) {
            new_edge_data[source - 1] = new_target_data_map;
        }
        else{
            new_edge_data[source] = new_target_data_map;
        }
    }
    this->_edge_data.clear();
    this->_edge_data = new_edge_data;
}

inline void UDataGraph::RelabelNode(const NodeId nodeId, const Label newLabel) {
    UGraphStruct::RelabelNode(nodeId, newLabel);
}


/// Get the node data for some node in the directed data _graph (by type)
/// \param node
/// \param type
/// \return
inline double UDataGraph::GetNodeData(NodeId node, const std::string &type) const {
    return GetNodeData(node, _node_data_names.at(type));
}

/// Get the node data for some node in the directed data _graph (by index)
/// \param node
/// \param index
/// \return
inline double UDataGraph::GetNodeData(NodeId node, int index) const {
    return _node_data[node][index];
}

/// Get the node data for some node in the directed data _graph (all data)
/// \param node
/// \return
inline const std::vector<double> &UDataGraph::GetNodeData(NodeId node) const {
    return _node_data[node];
}

inline bool UDataGraph::dijkstra(const UDataGraph &graph, NodeId src, double (*weight_function)(const UDataGraph&, EDGE), std::vector<double>& distances) {
    std::priority_queue p_queue = std::priority_queue<std::pair<double, NodeId>, std::vector<std::pair<double, NodeId>>, std::greater<>>();
    //initialize distances
    for (NodeId i = 0; i < graph.nodes(); i++)
        distances[i] = std::numeric_limits<double>::max();
    distances[src] = 0;
    p_queue.emplace(0.0, src);
    NodeId nodesFound = 0;
    // Find the shortest paths
    while(!p_queue.empty()){
        const std::pair<double, NodeId>& distance_node_pair = p_queue.top();
        ++nodesFound;
        NodeId nodeId = distance_node_pair.second;
        double distance = distance_node_pair.first;
        p_queue.pop();
        for (NodeId x : graph.get_neighbors(nodeId)) {
            if (distances[x] == std::numeric_limits<double>::max()) {
                distances[x] = distance + weight_function(graph, EDGE(nodeId, x));
                p_queue.emplace(distances[x], x);
            }
            else{
                distances[x] = std::min(distances[x], distance + weight_function(graph, EDGE(nodeId, x)));
            }
        }
    }
    return nodesFound == graph.nodes();
}

inline bool UDataGraph::dijkstra(const UDataGraph &graph, NodeId src, NodeId dest,
                                 double (*weight_function)(const UDataGraph &, EDGE), std::vector<NodeId>& path, std::vector<double>& distances, double& length) {
    length = std::numeric_limits<double>::max();
    std::priority_queue p_queue = std::priority_queue<std::pair<double, NodeId>, std::vector<std::pair<double, NodeId>>, std::greater<>>();
    std::unordered_map<NodeId, std::pair<NodeId, int>> predecessors;
    //initialize distances
    for (NodeId i = 0; i < graph.nodes(); i++)
        distances[i] = std::numeric_limits<double>::max();
    distances[src] = 0;
    predecessors[0] = std::make_pair(0, 1);
    p_queue.emplace(0.0, src);
    // Find the shortest paths
    bool found = false;
    while(!p_queue.empty()){
        const std::pair<double, NodeId>& distance_node_pair = p_queue.top();
        NodeId nodeId = distance_node_pair.second;
        double distance = distance_node_pair.first;
        if (nodeId == dest){
            found = true;
            length = 0;
            break;
        }
        p_queue.pop();
        for (NodeId x : graph.get_neighbors(nodeId)) {
            if (distances[x] == std::numeric_limits<double>::max()) {
                distances[x] = distance + weight_function(graph, EDGE(nodeId, x));
                p_queue.emplace(distances[x], x);
                predecessors[x] = std::make_pair(nodeId, predecessors[nodeId].second + 1);
            }
            else{
                if (distances[x] > distance + weight_function(graph, EDGE(nodeId, x))){
                    distances[x] = distance + weight_function(graph, EDGE(nodeId, x));
                    predecessors[x] = std::make_pair(nodeId, predecessors[nodeId].second + 1);
                }
            }
        }
    }
    if (found){
        NodeId currentNode = dest;
        path.resize(predecessors[dest].second);
        path[0] = src;
        int index = predecessors[dest].second - 1;
        while (currentNode != src){
            path[index] = currentNode;
            length += weight_function(graph, EDGE(predecessors[currentNode].first, currentNode));
            currentNode = predecessors[currentNode].first;
            --index;
        }

    }
    return found;
}

inline void UDataGraph::Init(const std::string &name, int size, int edges, int nodeFeatures, int edgeFeatures,
                             const std::vector<std::string> &nodeFeatureNames,
                             const std::vector<std::string> &edgeFeatureNames) {
    UGraphStruct::Init(name, size, edges, nodeFeatures, edgeFeatures, nodeFeatureNames, edgeFeatureNames);
    this->_node_data_size = (int) nodeFeatureNames.size();
    this->_edge_data_size = (int) edgeFeatureNames.size();

    for (int i = 0; i < this->_node_data_size; ++i) {
        this->_node_data_names.insert({nodeFeatureNames[i], i});
    }
    for (int i = 0; i < this->_edge_data_size; ++i) {
        this->_edge_data_names.insert({edgeFeatureNames[i], i});
    }
    this->_node_data.resize(nodes());
}

void UDataGraph::ReadNodeFeatures(double value, INDEX pos, const std::string &nodeFeatureName) {
    this->_node_data[pos].emplace_back(value);
    if (nodeFeatureName == "label") {
        this->_labels.emplace_back((int) value);
    }
}


inline void UDataGraph::WriteGraph(std::ofstream& Out, const SaveParams& saveParams){
    const std::string graphName = this->_name;
    unsigned int stringLength = graphName.length();
    GraphType Type = this->graphType;
    auto Size = this->nodes();
    unsigned int numFeatures = 0;
    std::vector<std::string> featureNames;

    Out.write((char *) (&stringLength), sizeof(stringLength));
    Out.write(graphName.c_str(), stringLength);
    Out.write((char *) (&Type), sizeof(GraphType));
    Out.write((char *) (&Size), sizeof(INDEX));

    numFeatures = this->_node_data_size;
    featureNames.resize(numFeatures);
    for (auto &[n, id] : this->_node_data_names) {
        featureNames[id] = n;
    }

    Out.write((char *) (&numFeatures), sizeof(unsigned int));
    for (int j = 0; j < numFeatures; ++j) {
        unsigned int nodeFeatureStringLength = featureNames[j].length();
        Out.write((char *) (&nodeFeatureStringLength), sizeof(nodeFeatureStringLength));
        Out.write(featureNames[j].c_str(), nodeFeatureStringLength);
    }

    numFeatures = 0;
    featureNames.clear();
    auto Edges = this->edges();
    Out.write((char *) (&Edges), sizeof(INDEX));

    numFeatures = this->_edge_data_size;
    featureNames.resize(numFeatures);
    for (auto &[e, id] : this->_edge_data_names) {
        featureNames[id] = e;
    }

    Out.write(reinterpret_cast<char *>(&numFeatures), sizeof(unsigned int));
    for (int j = 0; j < numFeatures; ++j) {
        unsigned int edgeFeatureStringLength = featureNames[j].length();
        Out.write(reinterpret_cast<char *>(&edgeFeatureStringLength), sizeof(edgeFeatureStringLength));
        Out.write(featureNames[j].c_str(), edgeFeatureStringLength);
    }
}

inline void UDataGraph::WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams){
    for (int j = 0; j < this->nodes(); ++j) {
        for (int k = 0; k < _node_data_size; ++k) {
            auto val = _node_data[j][k];
            Out.write(reinterpret_cast<char *>(&val), sizeof(double));
        }
    }
}

inline void UDataGraph::WriteEdges(std::ofstream& Out,const SaveParams& saveParams){
    INDEX Src = 0;
    for (auto const &edges: this->_graph) {
        for (auto Dst: edges) {
            if (Dst >= Src) {
                //std::cout << "Saving Edge " << Src << " " << Dst << std::endl;
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
                for (int j = 0; j < _edge_data_size; ++j) {
                    auto val = _edge_data[Src][Dst][j];
                    Out.write(reinterpret_cast<char *>(&val), sizeof(val));
                }
            }
        }
        ++Src;
    }
}

inline INDEX UDataGraph::AddNodes(const INDEX number) {
    return UGraphStruct::AddNodes(number);
}

inline INDEX UDataGraph::AddNodes(const INDEX number, const std::vector<Label> &labels) {
    return UGraphStruct::AddNodes(number, labels);
}

#endif //TESTGRAPHLIB_GRAPHLABELEDBASE_H
