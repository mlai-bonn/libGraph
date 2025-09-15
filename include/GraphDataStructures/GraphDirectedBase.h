//
// Created by florian on 23.10.23.
//

#ifndef GRAPH_DIRECTED_BASE_H
#define GRAPH_DIRECTED_BASE_H

#include "GraphBase.h"

struct DGraphStruct : public GraphStruct {
    //Constructors
    DGraphStruct()= default;
    explicit DGraphStruct(const std::string & graphPath, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "");
    explicit DGraphStruct(GraphStruct& graph);
    explicit DGraphStruct(const std::string& name, INDEX size);
    explicit DGraphStruct(const std::string&name, INDEX size, const Labels& labels);



    //Load and save
    void Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format = "");
    void Save(const SaveParams& saveParams) override;

    void Init(const std::string& name, int size, int edges, int nodeFeatures, int edgeFeatures, const std::vector<std::string>& nodeFeatureNames, const std::vector<std::string>& edgeFeatureNames) override;


    // Graph manipulation (functions that change the underlying graph data)
    // change behaviour of adding edges because of directedness
    INDEX AddNodes(INDEX number) override;
    INDEX AddNodes(INDEX number, const std::vector<Label>& labels) override;
    INDEX AddNodes(INDEX number, const std::vector<Label>& labels, const std::vector<std::vector<double>>& nodeData) override;
    // Edges
    bool AddEdge(NodeId source, NodeId destination) override;
    bool AddEdge(const std::pair<NodeId, NodeId>& edge) override;

    bool AddEdge(NodeId source, NodeId destination, bool check_existence) override;
    bool AddEdge(const std::pair<NodeId, NodeId>& edge, bool check_existence) override;

    bool RemoveEdge(NodeId source, NodeId destination) override;
    bool RemoveEdge(const std::pair<NodeId, NodeId>& edge) override;


    void ReadNodeFeatures(double value, INDEX pos, const std::string& nodeFeatureName) override;
    void WriteEdges(std::ofstream& Out,const SaveParams& saveParams) override;


    bool IsEdge(NodeId source, NodeId destination) const override;
    bool edge(NodeId source, NodeId destination, bool directed) const;
    static DGraphStruct GetBFSTree(const GraphStruct& graph, NodeId rootNodeId);
    INDEX out_degree(NodeId node);
    INDEX in_degree(NodeId node);

    bool operator==(const DGraphStruct &rhs) const {
        if ((GraphStruct) *this != (GraphStruct) rhs){
            return false;
        }
        if (_in_degrees != rhs._in_degrees){
            return false;
        }
        if (_out_degrees != rhs._out_degrees){
            return false;
        }
        return true;
    }

    //TODO out and in degree
    std::vector<NodeId> _in_degrees;
    std::vector<NodeId> _out_degrees;
};


/// Construct a directed _graph from an undirected _graph
/// \param graph
inline DGraphStruct::DGraphStruct(GraphStruct & graph) : GraphStruct(graph){
    for (INDEX i = 0; i < graph.nodes(); ++i) {
        this->_in_degrees.emplace_back(graph.degree(i));
        this->_out_degrees.emplace_back(graph.degree(i));
    }
}

/// Construct a directed _graph by certain size
/// \param name
/// \param size
inline DGraphStruct::DGraphStruct(const std::string& name, const INDEX size) : GraphStruct(name, size, Labels()){
    this->_in_degrees.resize(size, 0);
    this->_out_degrees.resize(size, 0);
}

/// Constructor of directed _graph with certain size and given labels
/// \param name
/// \param size
/// \param labels
inline DGraphStruct::DGraphStruct(const std::string& name, const INDEX size, const Labels &labels) : GraphStruct(name, size, labels) {
    this->_in_degrees.resize(size, 0);
    this->_out_degrees.resize(size, 0);
}

/// Read directed _graph from file
/// \param graphPath
/// \param relabeling
/// \param withLabels
/// \param labelPath
inline DGraphStruct::DGraphStruct(const std::string & graphPath, bool relabeling, bool withLabels, const std::string& labelPath){
    GraphStruct::Load(graphPath, relabeling, withLabels, labelPath);
}

inline INDEX DGraphStruct::out_degree(NodeId node) {
    return this->_out_degrees[node];
}



/// Ad an edge to a directed _graph
/// \param source
/// \param destination
/// \param check_existence
/// \return
inline bool DGraphStruct::AddEdge(const NodeId source, NodeId destination, const bool check_existence) {
    if (!check_existence || !IsEdge(source, destination)){
        this->_graph[source].emplace_back(destination);
        NodeId ElementId = (NodeId) this->_graph[source].size() - 1;
        while (ElementId > 0 && this->_graph[source][ElementId] < this->_graph[source][ElementId - 1]){
            std::swap(this->_graph[source][ElementId], this->_graph[source][ElementId - 1]);
            --ElementId;
        }
        ++this->_edges;
        ++this->_in_degrees[destination];
        ++this->_out_degrees[source];
        ++this->_degrees[destination];
        ++this->_degrees[source];
        this->maxDegree = std::max(this->maxDegree, static_cast<int>(std::max(this->degree(source), this->degree(destination))));
        return true;
    }
    return false;
}

inline bool DGraphStruct::AddEdge(const std::pair<NodeId, NodeId> &edge, const bool check_existence) {
    return DGraphStruct::AddEdge(edge.first, edge.second, check_existence);
}

inline bool DGraphStruct::RemoveEdge(const NodeId source, const NodeId destination) {
    if (const auto it_source = std::ranges::find(this->_graph[source], destination); it_source != this->_graph[source].end()) {
        this->_graph[source].erase(it_source);
        --_degrees[source];
        --_degrees[destination];
        --this->_in_degrees[destination];
        --this->_out_degrees[source];
        --this->_edges;
        return true;
    }
    return false;
}

inline bool DGraphStruct::RemoveEdge(const std::pair<NodeId, NodeId> &edge) {
    return DGraphStruct::RemoveEdge(edge.first, edge.second);
}

/// Check if edge in directed _graph exists
/// \param source
/// \param destination
/// \return
inline bool DGraphStruct::IsEdge(NodeId source, NodeId destination) const {
    return DGraphStruct::edge(source, destination, true);
}

/// Check if edge in directed _graph exists
/// \param source
/// \param destination
/// \param directed
/// \return
inline bool DGraphStruct::edge(NodeId source, NodeId destination, bool directed) const {
    if (directed){
        return std::binary_search(this->_graph[source].begin(), this->_graph[source].end(), destination);
    }
    else{
        if (std::binary_search(this->_graph[source].begin(), this->_graph[source].end(), destination)){
            return true;
        }
        return std::binary_search(this->_graph[destination].begin(), this->_graph[destination].end(), source);
    }
}

inline INDEX DGraphStruct::in_degree(NodeId node) {
    return _in_degrees[node];
}
///
/// \param graph
/// \param rootNodeId
/// \return
inline DGraphStruct DGraphStruct::GetBFSTree(const GraphStruct &graph, const NodeId rootNodeId) {
    DGraphStruct BFSTree = DGraphStruct(graph.GetName() + "BFS_Tree", graph.nodes(), graph.labels());
    std::vector<bool> VisitedNodes = std::vector<bool>(graph.nodes(), false);
    VisitedNodes[rootNodeId] = true;
    std::vector<NodeId> CurrentNodes;
    CurrentNodes.push_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        for (NodeId neighborId : graph.get_neighbors(NextNodeId)) {
            if (!VisitedNodes[neighborId]){
                VisitedNodes[neighborId] = true;
                CurrentNodes.insert(CurrentNodes.begin(), neighborId);
                BFSTree.AddEdge(NextNodeId, neighborId);
            }
        }
    }
    return BFSTree;
}

inline void DGraphStruct::Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format) {
    if (std::filesystem::is_regular_file(graphPath)) {
        int Version = 1;
        int graphId = 0;
        std::string extension = std::filesystem::path(graphPath).extension().string();
        if (extension == ".bgf" || extension == ".bgfs"){
            GraphStruct::Load(graphPath,relabeling,withLabels,labelPath);
        }
        else if (extension == ".bin") {
            std::ifstream In(graphPath, std::ios::in | std::ios::binary);
            std::string Name;
            unsigned int stringLength;
            In.read((char *) (&stringLength), sizeof(stringLength));
            Name.resize(stringLength);
            In.read((char *) Name.c_str(), stringLength);
            this->_name = Name;
            In.read((char *) (&this->graphType), sizeof(GraphType));
            In.read((char *) (&this->_nodes), sizeof(INDEX));
            In.read((char *) (&this->_edges), sizeof(INDEX));
            this->_in_degrees.resize(this->nodes());
            this->_out_degrees.resize(this->nodes());
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
                    ++this->_out_degrees[Src];
                    ++this->_in_degrees[Dst];
                }
                else{
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
            std::string graph_name = std::filesystem::path(graphPath).stem().string();
            this->_name = graph_name;
            NodeId src;
            NodeId dest;
            std::string a, b;
            std::string line;
            std::ifstream infile(graphPath);
            std::vector<std::pair<INDEX, INDEX>> graphEdges;
            std::set<INDEX> graphNodeIds;
            std::unordered_map<INDEX,INDEX> originalIdsToNodeIds;
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
            this->_degrees.resize(_nodes);
            this->_graph.resize(_nodes);
            this->_in_degrees.resize(_nodes);
            this->_out_degrees.resize(_nodes);
            INDEX nodeCounter = 0;
            for (auto x : graphNodeIds) {
                originalIdsToNodeIds.insert({x,nodeCounter});
                ++nodeCounter;
            }
            for (auto edge : graphEdges) {
                GraphStruct::AddEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second], true);
            }
            if (!labelPath.empty() && withLabels) {
                std::string label_extension = std::filesystem::path(labelPath).extension().string();
                std::ifstream label_file(labelPath);
                if (label_extension == ".vertexids"){
                    std::string label_line;
                    _labels = Labels(nodes(), 0);
                    Label label = 0;
                    while(std::getline(label_file, label_line)){
                        std::istringstream iss(label_line);
                        std::string token;
                        while (std::getline(iss, token, ' ')){
                            _labels[originalIdsToNodeIds[std::stoi(token)]] = label;
                        }
                        ++label;
                    }
                }
                else {
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
    if (CheckTree()) {
        this->graphType = GraphType::TREE;
    }
}

inline void DGraphStruct::Save(const SaveParams& saveParams) {
    std::string saveName;
    if (saveParams.graphPath.empty()) {
        if (!this->_path.empty()) {
            saveName = this->_path + this->_name;
        } else {
            return;
        }
    } else {
        saveName = saveParams.graphPath + this->_name;
    }
    if (!saveParams.Name.empty()) {
        saveName = saveParams.graphPath + saveParams.Name;
    }
    switch (saveParams.Format) {
        case GraphFormat::BGF: {
            GraphStruct::Save(saveParams);
            break;
        }
        case GraphFormat::BGFS: {
            GraphStruct::Save(saveParams);
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
                    Out.write((char *) (&Src), sizeof(INDEX));
                    Out.write((char *) (&Dst), sizeof(INDEX));
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
                    if (!FirstLine) {
                        Out << "\n";
                    }
                    Out << Src << " " << Dst;
                    FirstLine = false;
                }
                ++Src;
            }
        }
        case GraphFormat::PEREGRINE_DATA:
            break;
        case GraphFormat::PEREGRINE_SMALL:
            break;
    }
}

inline void DGraphStruct::Init(const std::string &name, int size, int edges, int nodeFeatures, int edgeFeatures,
                        const std::vector<std::string> &nodeFeatureNames,
                        const std::vector<std::string> &edgeFeatureNames) {
    GraphStruct::Init(name, size, edges, nodeFeatures, edgeFeatures, nodeFeatureNames, edgeFeatureNames);
    //Create _graph
    this->_in_degrees.resize(_nodes);
    this->_out_degrees.resize(_nodes);
}

inline INDEX DGraphStruct::AddNodes(const INDEX number, const std::vector<Label>& labels) {
    // Emplace back in and out degrees
    this->_in_degrees.emplace_back(0);
    this->_out_degrees.emplace_back(0);
    return GraphStruct::AddNodes(number, labels);
}

inline INDEX DGraphStruct::AddNodes(const INDEX number, const std::vector<Label> &labels,
    const std::vector<std::vector<double>> &nodeData) {
    return GraphStruct::AddNodes(number, labels, nodeData);
}

inline INDEX DGraphStruct::AddNodes(const INDEX number) {
    // Emplace back in and out degrees
    this->_in_degrees.emplace_back(0);
    this->_out_degrees.emplace_back(0);
    return GraphStruct::AddNodes(number);
}

inline bool DGraphStruct::AddEdge(const NodeId source, const NodeId destination) {
    return DGraphStruct::AddEdge(source, destination, true);
}

inline bool DGraphStruct::AddEdge(const std::pair<NodeId, NodeId> &edge) {
    return DGraphStruct::AddEdge(edge.first, edge.second);
}

inline void DGraphStruct::ReadNodeFeatures(double value, INDEX pos, const std::string &nodeFeatureName) {
    GraphStruct::ReadNodeFeatures(value, pos, nodeFeatureName);
}

inline void DGraphStruct::WriteEdges(std::ofstream& Out, const SaveParams& saveParams){
    INDEX Src = 0;
    for (auto const &edges: this->_graph) {
        for (auto Dst: edges) {
            if (saveParams.Format == GraphFormat::BGF) {
                Out.write((char *) (&Src), sizeof(INDEX));
                Out.write((char *) (&Dst), sizeof(INDEX));
            }
            else{
                auto int_Src = (unsigned int) Src;
                Out.write((char *) (&int_Src), sizeof(unsigned int));
                auto int_Dst = (unsigned int) Dst;
                Out.write((char *) (&int_Dst), sizeof(unsigned int));
            }
            unsigned int edgeFeatures = 0;
            for (int j = 0; j < edgeFeatures; ++j) {
                //Out.write((char *) (&Src), sizeof(INDEX));
            }
        }
        ++Src;
    }
}




#endif //GRAPH_DIRECTED_BASE_H
