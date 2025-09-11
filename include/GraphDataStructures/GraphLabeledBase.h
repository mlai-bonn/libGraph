//
// Created by florian on 23.10.23.
//

#ifndef TESTGRAPHLIB_GRAPHLABELEDBASE_H
#define TESTGRAPHLIB_GRAPHLABELEDBASE_H
#include "GraphDirectedBase.h"

struct DDataGraph : public DGraphStruct{
public:
    DDataGraph();
    explicit DDataGraph(const std::string & graphPath, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "");
    void Load(const std::string & graphPath, bool relabeling, bool withLabels, const std::string& labelPath, const std::string& format = "");
    void Save(const SaveParams& saveParams) override;

    void Init(const std::string& name, int size, int edges, int nodeFeatures, int edgeFeatures, const std::vector<std::string>& nodeFeatureNames, const std::vector<std::string>& edgeFeatureNames) override;
    void ReadNodeFeatures(double value, INDEX pos, const std::string& nodeFeatureName) override;
    bool ReadEdges(INDEX Src, INDEX Dst, std::vector<double>& edgeData) override;

    void WriteGraph(std::ofstream& Out, const SaveParams& saveParams) override;
    void WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams) override;
    void WriteEdges(std::ofstream& Out,const SaveParams& saveParams) override;

    NodeId add_node(INDEX number = 1, const Labels* labels = nullptr) override;


    double get_edge_data(const EDGE & edge, const std::string& type) const;
    double get_edge_data(const EDGE & edge, int index) const;
    const std::vector<double>& get_edge_data(const EDGE& edge) const;
    double get_node_data(NodeId node, const std::string& type) const;
    double get_node_data(NodeId node, int index) const;
    const std::vector<double>& get_node_data(NodeId node) const;
    void add_edge_data(const EDGE & edge, const std::string& type, double data);
    void add_edge_data(const EDGE & edge, int index, double data);
    void add_edge_data(const EDGE & edge, std::vector<double>& data);
    bool add_edge(NodeId source, NodeId destination, std::vector<double>& data);

    // set node data names
    void set_node_data_names(const std::unordered_map<std::string, int>& nodeDataNames);
    // set edge data names
    void set_edge_data_names(const std::unordered_map<std::string, int>& edgeDataNames);


    static bool dijkstra(const DDataGraph& graph, NodeId src, double (*weight_function)(const DDataGraph&, EDGE), std::vector<double>& distances);
    static bool dijkstra(const DDataGraph& graph, NodeId src, NodeId dest, double (*weight_function)(const DDataGraph&, EDGE), std::vector<NodeId>& path, std::vector<double>& distances, double& length);

    bool operator==(const DDataGraph &rhs) const {
        if ((DGraphStruct) *this != (DGraphStruct) rhs){
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
inline DDataGraph::DDataGraph() = default;

/// Read directed data _graph from file TODO binary extension
/// \param graphPath path of the _graph file
/// \param relabeling if true node labels start with 0 and end with size-1, otherwise use original labels TODO check this
/// \param withLabels if considering node attributes
/// \param labelPath path to node data
inline DDataGraph::DDataGraph(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath){
    DDataGraph::Load(graphPath, relabeling, withLabels,labelPath);
}

inline void DDataGraph::Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format) {
    int Version = 1;
    int graphId = 0;
    if (std::filesystem::is_regular_file(graphPath)) {
        std::string extension = std::filesystem::path(graphPath).extension().string();
        if (extension == ".bgf"){
            DGraphStruct::Load(graphPath, relabeling,withLabels,labelPath);
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
            this->_in_degrees.resize(_nodes);
            this->_out_degrees.resize(_nodes);
            INDEX nodeCounter = 0;
            for (auto x: graphNodeIds) {
                originalIdsToNodeIds.insert({x, nodeCounter});
                this->IdsToOriginalIds.insert({nodeCounter, x});
                ++nodeCounter;
            }
            for (auto edge: graphEdges) {
                DDataGraph::add_edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second],
                                     original_edge_data[edge.first][edge.second]);
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
                    labelMap = GraphFunctions::GetGraphLabelMap(_labels);
                    labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
                    //TODO why this line this->numLabels = (this->numLabels == -1) ? static_cast<int>(labelMap.size()) : this->numLabels;
                    this->_numLabels = static_cast<int>(labelMap.size());
                    if (this->_numLabels >= 10) {
                        labelType = LABEL_TYPE::LABELED_DENSE;
                    } else {
                        labelType = LABEL_TYPE::LABELED_SPARSE;
                    }
                    UpdateGraphLabels(labelType);
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

inline void DDataGraph::Save(const SaveParams& saveParams) {
    DGraphStruct::Save(saveParams);
}

/// Get the data of an edge using the data type
/// \param edge
/// \param type
/// \return
inline double DDataGraph::get_edge_data(const EDGE & edge, const std::string& type) const {
    int index = _edge_data_names.at(type);
    return _edge_data.at(edge.first).at(edge.second)[index];
}

/// Get the data of an edge using the data index
/// \param edge
/// \param index
/// \return
inline double DDataGraph::get_edge_data(const EDGE &edge, int index) const {
    return _edge_data.at(edge.first).at(edge.second)[index];
}

/// Get edge data for some edge in a directed data _graph
/// \param edge
/// \return
inline const std::vector<double>& DDataGraph::get_edge_data(const EDGE & edge) const {
    return _edge_data.at(edge.first).at(edge.second);
}


/// Add the data of an edge using the type
/// \param edge
/// \param type
/// \param data
inline void DDataGraph::add_edge_data(const EDGE &edge, const std::string &type, double data) {
    add_edge_data(edge, _edge_data_names[type], data);
}

/// Add the data of an edge using the edge index
/// \param edge
/// \param index
/// \param data
inline void DDataGraph::add_edge_data(const EDGE &edge, int index, double data) {
    _edge_data[edge.first][edge.second].resize(this->_edge_data_size);
    _edge_data[edge.first][edge.second][index] = data;
}

/// Add the data of an edge
/// \param edge
/// \param data
inline void DDataGraph::add_edge_data(const EDGE &edge, std::vector<double> &data) {
    _edge_data[edge.first][edge.second] = data;
}

/// Add an data edge to a directed data _graph
/// \param source
/// \param destination
/// \param data
/// \return
inline bool DDataGraph::add_edge(NodeId source, NodeId destination, std::vector<double>& data) {
    add_edge_data(EDGE{source, destination}, data);
    return DGraphStruct::add_edge(source, destination);
}

inline void DDataGraph::set_node_data_names(const std::unordered_map<std::string, int> &nodeDataNames) {
    this->_node_data_names = nodeDataNames;
    this->_node_data_size = (int) nodeDataNames.size();
}

inline void DDataGraph::set_edge_data_names(const std::unordered_map<std::string, int> &edgeDataNames) {
    this->_edge_data_names = edgeDataNames;
    this->_edge_data_size = (int) edgeDataNames.size();
}

inline NodeId DDataGraph::add_node(INDEX number, const Labels *labels) {
    const NodeId new_node_id =  DGraphStruct::add_node(number, labels);
    this->_in_degrees.emplace_back(0);
    this->_out_degrees.emplace_back(0);
    return new_node_id;
}


/// Get the node data for some node in the directed data _graph (by type)
/// \param node
/// \param type
/// \return
inline double DDataGraph::get_node_data(NodeId node, const std::string &type) const {
    return get_node_data(node, _node_data_names.at(type));
}

/// Get the node data for some node in the directed data _graph (by index)
/// \param node
/// \param index
/// \return
inline double DDataGraph::get_node_data(NodeId node, int index) const {
    return _node_data[node][index];
}

/// Get the node data for some node in the directed data _graph (all data)
/// \param node
/// \return
inline const std::vector<double> &DDataGraph::get_node_data(NodeId node) const {
    return _node_data[node];
}

inline bool DDataGraph::dijkstra(const DDataGraph &graph, NodeId src, double (*weight_function)(const DDataGraph&, EDGE), std::vector<double>& distances) {
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

inline bool DDataGraph::dijkstra(const DDataGraph &graph, NodeId src, NodeId dest,
                                 double (*weight_function)(const DDataGraph &, EDGE), std::vector<NodeId>& path, std::vector<double>& distances, double& length) {
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

inline void DDataGraph::Init(const std::string &name, int size, int edges, int nodeFeatures, int edgeFeatures,
                             const std::vector<std::string> &nodeFeatureNames,
                             const std::vector<std::string> &edgeFeatureNames) {
    DGraphStruct::Init(name, size, edges, nodeFeatures, edgeFeatures, nodeFeatureNames, edgeFeatureNames);
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

void DDataGraph::ReadNodeFeatures(double value, INDEX pos, const std::string &nodeFeatureName) {
    this->_node_data[pos].emplace_back(value);
    if (nodeFeatureName == "label") {
        this->_labels.emplace_back((int) value);
    }
}

bool DDataGraph::ReadEdges(INDEX Src, INDEX Dst, std::vector<double> &edgeData) {
    return this->add_edge(Src, Dst, edgeData);
}


inline void DDataGraph::WriteGraph(std::ofstream& Out, const SaveParams& saveParams){
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

inline void DDataGraph::WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams){
    for (int j = 0; j < this->nodes(); ++j) {
        for (int k = 0; k < _node_data_size; ++k) {
            auto val = _node_data[j][k];
            Out.write(reinterpret_cast<char *>(&val), sizeof(double));
        }
    }
}

inline void DDataGraph::WriteEdges(std::ofstream& Out,const SaveParams& saveParams){
    INDEX Src = 0;
    for (auto const &edges: this->_graph) {
        for (auto Dst: edges) {
            Out.write(reinterpret_cast<char *>(&Src), sizeof(INDEX));
            Out.write(reinterpret_cast<char *>(&Dst), sizeof(INDEX));
            for (int j = 0; j < _edge_data_size; ++j) {
                auto val = _edge_data[Src][Dst][j];
                Out.write(reinterpret_cast<char *>(&val), sizeof(val));
            }
        }
        ++Src;
    }
}

#endif //TESTGRAPHLIB_GRAPHLABELEDBASE_H
