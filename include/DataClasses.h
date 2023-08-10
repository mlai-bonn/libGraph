//
// Created by Florian on 15.04.2021.
//

#ifndef HOPS_DATACLASSES_H
#define HOPS_DATACLASSES_H
#include <vector>
#include "typedefs.h"
#include "GraphFunctions.h"
#include <set>
#include <random>
#include <unordered_map>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>

enum class GraphType{
    GENERAL,
    TREE,
    OUTERPLANAR,
};

enum class GraphFormat{
    BGF,
    BGFS,
    BINARY,
    EDGES,
    PEREGRINE_DATA,
    PEREGRINE_SMALL,
    DIMACS,
    AIDS,
};

struct DGraphStruct;

struct SaveParams{
    std::string graphPath = "";
    std::string Name = "";
    GraphFormat Format = GraphFormat::BGFS;
    bool Labeled = false;
};

struct GraphStruct{
public:
    GraphStruct()= default;
    explicit GraphStruct(const std::string & graphPath, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "", const std::string& formate = "", const std::string& search_name = "");
    GraphStruct(NodeId size, const Labels& labels);

    void Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format = "", const std::string& search_name = "");
    virtual void Save(const SaveParams& saveParams);

    virtual void Init(const std::string& name, int size, int edges, int nodeFeatures, int edgeFeatures, const std::vector<std::string>& nodeFeatureNames, const std::vector<std::string>& edgeFeatureNames);
    virtual void ReadNodeFeatures(double value, int pos, const std::string& nodeFeatureName);
    virtual bool ReadEdges(INDEX Src, INDEX Dst, std::vector<double>& edgeData);

    virtual void WriteGraph(std::ofstream& Out, const SaveParams& saveParams);
    virtual void WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams);
    virtual void WriteEdges(std::ofstream& Out,const SaveParams& saveParams);

    void SetName(const std::string& name);
    std::string GetName() const{return _name;};
    void SetNumLabels(int numLabels){_numLabels = numLabels;};
    int GetNumLabels() const{return _numLabels;};
    std::string GetPath(){return _path;};

    void write_graph_nodes(const std::string& graphPath, const std::string& fileName, const Nodes& nodes);
    void CreateGraph(NodeId size, const Labels& labels);
    void UpdateGraphLabels(LABEL_TYPE label);

    const std::vector<Nodes>& graph() const {return _graph;};
    INDEX nodes() const;
    INDEX edges() const;
    void set_graph(const std::vector<Nodes>& nodes);
    void set_edges(INDEX edges){_edges=edges;};
    INDEX degree(NodeId nodeId) const;
    INDEX degree(NodeId nodeId, Label label);
    const std::vector<INDEX>& degreeByLabel(NodeId nodeId);
    const std::vector<Label>& labels() const;
    Label label(NodeId nodeId) const;
    void set_labels(const Labels *Labels);
    const Nodes& get_neighbors(NodeId nodeId) const;
    Nodes& neighbors(NodeId nodeId);
    NodeId neighbor(NodeId nodeId, INDEX neighborIdx) const;
    NodeId neighbor(NodeId nodeId, INDEX neighborIdx, Label label);
    NodeId random_neighbor_in_range(NodeId nodeId, INDEX minIdx, std::mt19937_64& gen);
    virtual bool edge(NodeId source, NodeId destination) const;
    virtual bool add_edge(NodeId source, NodeId destination);
    NodeId add_node(INDEX number = 1, Labels* labels = nullptr);
    //Get get_neighbors by []
    const Nodes& operator[](NodeId nodeId){return _graph[nodeId];};
    void sortNeighborIds();
    bool has_neighbor_label(NodeId nodeId, Label label);

    bool comp_degree(NodeId i, NodeId j) const{
        return this->degree(i) < this->degree(j);
    };

    void degree_sort(Nodes& nodes){
        std::sort(nodes.begin(), nodes.end(), [this](NodeId l,NodeId r){return comp_degree(l,r);});
    }

    static bool ReadBGF(const std::string& extension, std::ifstream& In, int& saveVersion, int& graphNumber, std::vector<std::string>& graphsNames,
                         std::vector<GraphType>& graphsTypes,
                         std::vector<INDEX>& graphsSizes,
                         std::vector<std::vector<std::string>>& graphsNodeFeatureNames,
                         std::vector<INDEX>& graphsEdges,
                         std::vector<std::vector<std::string>>& graphsEdgeFeatureNames);
    void InitLabels();
    template<typename T>
    static void LoadGraphsFromPath(const std::string& graphPath, const std::string& labelPath, std::vector<T> &graphs, const std::string & searchName = "", const std::string & extension = "", bool sort = true, std::set<int>* graphSizes = nullptr, int patternNum = -1);
    static void Convert(const std::string & graphPath, const std::string& Name = "", GraphFormat Format = GraphFormat::BGF, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "", bool Labeled = false);

    /// Converts all graphs in a directory using the same extension to the given format
    /// \param path
    /// \param Format
    static void Convert(const std::string & path, GraphFormat Format = GraphFormat::BGF, const std::string & extension = ".txt");
    bool IsTree() const;
    bool CheckSpanningTree(const GraphStruct& spanningTree) const;

    static void BFSDistances(const GraphStruct &graph, INDEX root, std::vector<INDEX> &distances);
    static void GetComponents(const GraphStruct& graph, std::vector<INDEX> & components);
    static void GetComponents(const GraphStruct& graph, std::vector<Nodes> & components);
    static void GetLargestComponent(const GraphStruct& graph, Nodes& nodes);
    static GraphStruct GetLargestComponent(const GraphStruct& graph);
    static GraphStruct SubGraph(const GraphStruct& graph, const Nodes& nodeIds);
    static bool IsConnected(const GraphStruct& graph);
    static bool ReachableNodes(const GraphStruct &graph, NodeId root, std::vector<INDEX> & reachability, INDEX Id, INDEX& number);

    static void DFS(const GraphStruct &graph, GraphStruct &tree, Nodes& nodeOrder, NodeId rootNodeId = -1, int seed = 0);
    static void OptOrdering(const GraphStruct &graph, Nodes& nodeOrder);
    static void ReorderGraph(GraphStruct &graph, const Nodes& nodeOrder);

    //Static functions
    //compare the given labeled degree vector with the labeled degree vector of the node with Id:nodeId of this graph
    bool isBigger(const std::vector<INDEX>& labeledDegree, NodeId nodeId);

    //Iterators
    struct NodeIterator{
        // Prefix increment
        void operator++() { ++nodeId;};
        const Nodes& operator*() const { return _graph->get_neighbors(nodeId); }
        friend bool operator== (const NodeIterator& a, INDEX b) { return a.nodeId == b;};
        friend bool operator!= (const NodeIterator& a, INDEX b) { return a.nodeId != b;};

        const GraphStruct* _graph{};
        NodeId nodeId{};
    };
    [[nodiscard]] NodeIterator begin() const { return NodeIterator{this, 0}; }
    [[nodiscard]] INDEX end() const   { return this->nodes(); }

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
        if (_numLabels != rhs._numLabels){
            return false;
        }
        if (isTree != rhs.isTree){
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
    bool isTree = false;

    //for labels
    std::unordered_map<Label, Nodes> labelMap{};
    std::unordered_map<Label, INDEX> labelFrequencyMap{};


protected:
    // General variables
    std::string _name;
    std::string _path;
    GraphType graphType = GraphType::GENERAL;
    INDEX _nodes = 0;
    INDEX _edges = 0;
    std::vector<Nodes> _graph;
    std::vector<INDEX> _degrees;

    //store original node Ids
    std::unordered_map<INDEX ,INDEX> IdsToOriginalIds;

    //for labels
    int _numLabels = -1;
    Labels _labels = Labels();
    std::vector<std::unordered_map<Label, INDEX>> labeledDegreeMap{};
    std::vector<std::vector<INDEX>> labeledDegreeVector{};
    std::vector<std::vector<Nodes>> labeledVectorGraph{};
    std::vector<std::unordered_map<Label, Nodes>> labeledMapGraph{};

public :
    [[deprecated]]
    void Save_0(const std::string & graphPath = "", GraphFormat Format = GraphFormat::BINARY, bool Labeled = false, bool OnlyGraph = true, const std::string& Name = "") const;


};

struct DGraphStruct : public GraphStruct{
    //Constructors
    DGraphStruct()= default;
    explicit DGraphStruct(const std::string & graphPath, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "");
    explicit DGraphStruct(GraphStruct& graph);
    explicit DGraphStruct(INDEX size);
    explicit DGraphStruct(INDEX size, const Labels& labels);

    //Load and save
    void Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format = "");
    void Save(const SaveParams& saveParams) override;

    void Init(const std::string& name, int size, int edges, int nodeFeatures, int edgeFeatures, const std::vector<std::string>& nodeFeatureNames, const std::vector<std::string>& edgeFeatureNames) override;
    void ReadNodeFeatures(double value, int pos, const std::string& nodeFeatureName) override;
    bool ReadEdges(INDEX Src, INDEX Dst, std::vector<double>& edgeData) override;

    void WriteEdges(std::ofstream& Out,const SaveParams& saveParams) override;

    bool edge(NodeId source, NodeId destination) const override;
    bool edge(NodeId source, NodeId destination, bool directed) const;
    bool add_edge(NodeId source, NodeId destination) override;
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


struct DDataGraph : public DGraphStruct{
public:
    DDataGraph();
    explicit DDataGraph(const std::string & graphPath, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "");
    void Load(const std::string & graphPath, bool relabeling, bool withLabels, const std::string& labelPath, const std::string& format = "");
    void Save(const SaveParams& saveParams) override;

    void Init(const std::string& name, int size, int edges, int nodeFeatures, int edgeFeatures, const std::vector<std::string>& nodeFeatureNames, const std::vector<std::string>& edgeFeatureNames) override;
    void ReadNodeFeatures(double value, int pos, const std::string& nodeFeatureName) override;
    bool ReadEdges(INDEX Src, INDEX Dst, std::vector<double>& edgeData) override;

    void WriteGraph(std::ofstream& Out, const SaveParams& saveParams) override;
    void WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams) override;
    void WriteEdges(std::ofstream& Out,const SaveParams& saveParams) override;


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

template <typename T>
class GraphData {
public:
    GraphData()= default;
    explicit GraphData(const std::string& graphPath, const std::string& labelPath = "", const std::string& searchName = "", const std::string& extension = "", bool sort = true, std::set<int>* graphSizes = nullptr, int patternNum = -1);
    explicit GraphData(GraphFormat graphFormat, const std::string& graphPath, const std::string& search_name = "", bool sort = true);

    void Load(const std::string &graphPath);
    void Load(std::vector<T>& graphs, const std::string &graphPath, GraphFormat graphFormat = GraphFormat::BGFS);

    void LoadBGF(std::vector<T>& graphs, const std::string &graphPath, GraphFormat graphFormat);
    void LoadAIDS(std::vector<T>& graphs, const std::string &graphPath);

    void Save(const SaveParams& saveParams = SaveParams());

    void add_graph(const std::string& graphPath, bool withLabels = false);
    void add(const T& graph);
    void add(const std::vector<T>& graphs);
    void add(const std::string& graphPath, const std::string& labelPath = "", const std::string& searchName = "", const std::string& extension = ".edges", bool sort = true, std::set<int>* graphSizes = nullptr, int patternNum = -1);

    std::vector<T> graphData;
    T& operator[](size_t index){return graphData[index];};
    [[nodiscard]] INDEX size() const{return static_cast<INDEX>(graphData.size());};
};

template<typename T>
void GraphData<T>::LoadBGF(std::vector<T> &graphs, const std::string &graphPath, GraphFormat graphFormat) {
    int saveVersion = 1;
    int graphNumber;
    std::vector<std::string> graphsNames;
    std::vector<GraphType> graphsTypes;
    std::vector<INDEX> graphsSizes;
    std::vector<std::vector<std::string>> graphsNodeFeatureNames;
    std::vector<INDEX> graphsEdges;
    std::vector<std::vector<std::string>> graphsEdgeFeatureNames;
    std::ifstream In(graphPath, std::ios::in | std::ios::binary);
    std::string extension = ".bgf";
    if (graphFormat == GraphFormat::BGFS){
        extension = ".bgfs";
    }

    if (GraphStruct::ReadBGF(extension, In, saveVersion, graphNumber, graphsNames, graphsTypes, graphsSizes,
                             graphsNodeFeatureNames, graphsEdges, graphsEdgeFeatureNames)) {
        for (int i = 0; i < graphNumber; ++i) {
            graphs.emplace_back();
            T& graph = graphs.back();
            //Create graph
            graph.Init(graphsNames[i], (int) graphsSizes[i], (int) graphsEdges[i],
                       (int) graphsNodeFeatureNames[i].size(), (int) graphsEdgeFeatureNames[i].size(),
                       graphsNodeFeatureNames[i], graphsEdgeFeatureNames[i]);
            //Read the nodes
            for (int j = 0; j < graphsSizes[i]; ++j) {
                for (int k = 0; k < graphsNodeFeatureNames[i].size(); ++k) {
                    double val = 0;
                    if (extension == ".bgf") {
                        In.read((char *) (&val), sizeof(double));
                    } else if (extension == ".bgfs") {
                        unsigned int int_val;
                        In.read((char *) (&int_val), sizeof(unsigned int));
                        val = int_val;
                    }
                    graph.ReadNodeFeatures(val, j, graphsNodeFeatureNames[i][k]);

                }
            }
            //Read the edges
            INDEX added_edges = 0;
            INDEX original_edges = graph.edges();
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
                if (graph.ReadEdges(Src, Dst, edgeData)) {
                    ++added_edges;
                }
            }
            graph.set_edges(added_edges);
            graph.InitLabels();
        }
    }
    In.close();

}

template<typename T>
void GraphData<T>::LoadAIDS(std::vector<T> &graphs, const std::string &graphPath) {

    NodeId src;
    NodeId dest;
    std::string a, b, c;
    std::string line;
    std::ifstream infile(graphPath);
    std::vector<std::pair<INDEX, INDEX>> graphEdges;
    Labels edgeLabels;
    std::set<INDEX> graphNodeIds;
    std::unordered_map<INDEX, INDEX> originalIdsToNodeIds;
    Labels nodeLabels;
    unsigned int num_nodes = 0;
    unsigned int num_edges = 0;
    bool label_line = false;
    bool edge_line = false;
    int graph_counter = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        if (edge_line){
            iss >> a;
            iss >> b;
            iss >> c;
            graphEdges.emplace_back(std::stoi(a), std::stoi(b));
            graphNodeIds.emplace(std::stoi(a));
            graphNodeIds.emplace(std::stoi(b));
            edgeLabels.emplace_back(std::stoi(c));
            ++num_edges;
            edge_line = false;

            graphs.emplace_back();
            T& graph = graphs.back();
            //Create graph
            graph.Init(std::to_string(graph_counter), num_nodes, 0,1, 1,{"node_label"}, {"edge_label"});

            INDEX nodeCounter = 0;
            for (auto x: graphNodeIds) {
                originalIdsToNodeIds.insert({x, nodeCounter});
                ++nodeCounter;
            }
            unsigned int num_edges_duplicates = 0;
            for (auto edge: graphEdges) {
                if (!graph.add_edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (graph.edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            if (num_edges_duplicates > 0) {
                std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base graph)!"
                          << std::endl;
                std::cout << graphEdges.size() - num_edges_duplicates << " edges are remaining" << std::endl;
            }

            int max_label = *std::max(nodeLabels.begin(), nodeLabels.end());
            graph.set_labels(&nodeLabels);
            graph.InitLabels();

            graphEdges.clear();
            edgeLabels.clear();
            graphNodeIds.clear();
            originalIdsToNodeIds.clear();
            nodeLabels.clear();
            num_nodes = 0;
            num_edges = 0;

            ++graph_counter;
            continue;
        }
        iss >> a;
        if (a == "#" || a == "$") {
            label_line = true;
            continue;
        }
        if (label_line){
            ++num_nodes;
            nodeLabels.emplace_back(std::stoi(a));
            while (iss >> a){
                ++num_nodes;
                nodeLabels.emplace_back(std::stoi(a));
            }
            label_line = false;
            edge_line = true;
            continue;
        }

        iss >> b;
        // Check if file is of format
        // num_nodes
        // node_idA node_idB
        // ...
            src = std::stoull(a);
            dest = std::stoull(b);
            graphEdges.emplace_back(src, dest);
            graphNodeIds.emplace(src);
            graphNodeIds.emplace(dest);
    }
//    graph.nodes() = graphNodeIds.size();
//    NodeId num_edges = graphEdges.size();
//    this->_degrees.resize(graph.nodes());
//    this->_graph.resize(graph.nodes());
//    INDEX nodeCounter = 0;
//    for (auto x: graphNodeIds) {
//        originalIdsToNodeIds.insert({x, nodeCounter});
//        ++nodeCounter;
//    }
//    unsigned int num_edges_duplicates = 0;
//    for (auto edge: graphEdges) {
//        if (!GraphStruct::add_edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
//            std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
//                      << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
//            if (this->edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
//                std::cout << "It already exists!" << std::endl;
//                ++num_edges_duplicates;
//            }
//        }
//    }
//    std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base graph)!"
//              << std::endl;
//    std::cout << graphEdges.size() - num_edges_duplicates << " edges are remaining" << std::endl;
//    if (!labelPath.empty() && withLabels) {
//        std::string label_extension = std::filesystem::path(labelPath).extension().string();
//        std::ifstream label_file(labelPath);
//        if (label_extension == ".vertexids") {
//            std::string label_line;
//            _labels = Labels(nodes(), 0);
//            Label label = 0;
//            while (std::getline(label_file, label_line)) {
//                std::istringstream iss(label_line);
//                std::string token;
//                while (std::getline(iss, token, ' ')) {
//                    _labels[originalIdsToNodeIds[std::stoi(token)]] = label;
//                }
//                ++label;
//            }
//        } else {
//            try {
//                NodeId id;
//                Label label;
//                while (label_file >> id >> label) {
//                    _labels.push_back(label);
//                }
//            }
//            catch (...) {
//            }
//        }
//        if (this->_labels.size() == this->_graph.size()) {
//            labelMap = GraphFunctions::GetGraphLabelMap(_labels);
//            labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
//            this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size()) : this->_numLabels;
//            if (this->_numLabels >= 10) {
//                labelType = LABEL_TYPE::LABELED_DENSE;
//            } else {
//                labelType = LABEL_TYPE::LABELED_SPARSE;
//            }
//            UpdateGraphLabels(labelType);
//        } else {
//            //TODO throw exception
//        }
//    }
}

template<typename T>
void GraphData<T>::Save(const SaveParams& saveParams) {
    std::string saveName;
    if (saveParams.graphPath.empty() && this->graphData.size() > 0) {
        saveName += this->graphData[0].GetPath();
    }
    else{
        saveName += saveParams.graphPath;
    }
    if (!saveParams.Name.empty()) {
        saveName += saveParams.Name;
    }
    else {
        saveName += "GraphData";
    }
    switch (saveParams.Format) {
        case GraphFormat::BGF: {
            int SaveVersion = 1;
            int graphNumber = (int) this->graphData.size();
            std::ofstream Out(saveName + ".bgf", std::ios::out | std::ios::binary);
            Out.write((char *) (&SaveVersion), sizeof(int));
            Out.write((char *) (&graphNumber), sizeof(int));
            for (auto &graph: this->graphData) {
                graph.WriteGraph(Out, saveParams);
            }

            for (auto &graph: this->graphData) {
                graph.WriteNodeFeatures(Out, saveParams);
                graph.WriteEdges(Out, saveParams);
            }
            Out.close();
            break;
        }
        case GraphFormat::BGFS: {
            int SaveVersion = 1;
            int graphNumber = (int) this->graphData.size();
            std::ofstream Out(saveName + ".bgfs", std::ios::out | std::ios::binary);
            Out.write((char *) (&SaveVersion), sizeof(int));
            Out.write((char *) (&graphNumber), sizeof(int));
            for (auto &graph: this->graphData) {
                graph.WriteGraph(Out, saveParams);
            }
            for (auto &graph: this->graphData) {
                graph.WriteNodeFeatures(Out, saveParams);
                graph.WriteEdges(Out, saveParams);
            }
            Out.close();
            break;
        }
        case GraphFormat::BINARY:
            break;
        case GraphFormat::EDGES:
            break;
        case GraphFormat::PEREGRINE_DATA:
            break;
        case GraphFormat::PEREGRINE_SMALL:
            break;
    }
}

template<typename T>
void GraphData<T>::Load(const std::string &graphPath) {
    if (std::filesystem::is_regular_file(graphPath)) {
        std::string extension = std::filesystem::path(graphPath).extension().string();
        if (extension == ".bgf" || extension == ".bgfs") {
            int saveVersion = 1;
            int graphNumber;
            std::vector<std::string> graphsNames;
            std::vector<GraphType> graphsTypes;
            std::vector<INDEX> graphsSizes;
            std::vector<std::vector<std::string>> graphsNodeFeatureNames;
            std::vector<INDEX> graphsEdges;
            std::vector<std::vector<std::string>> graphsEdgeFeatureNames;
            std::ifstream In(graphPath, std::ios::in | std::ios::binary);

            if (GraphStruct::ReadBGF(extension, In, saveVersion, graphNumber, graphsNames, graphsTypes, graphsSizes,
                                     graphsNodeFeatureNames, graphsEdges, graphsEdgeFeatureNames)) {
                for (int i = 0; i < graphNumber; ++i) {
                    this->graphData.emplace_back();
                    T& graph = this->graphData.back();
                    //Create graph
                    graph.Init(graphsNames[i], (int) graphsSizes[i], (int) graphsEdges[i],
                                   (int) graphsNodeFeatureNames[i].size(), (int) graphsEdgeFeatureNames[i].size(),
                                   graphsNodeFeatureNames[i], graphsEdgeFeatureNames[i]);
                    //Read the nodes
                    for (int j = 0; j < graphsSizes[i]; ++j) {
                        for (int k = 0; k < graphsNodeFeatureNames[i].size(); ++k) {
                            double val = 0;
                            if (extension == ".bgf") {
                                In.read((char *) (&val), sizeof(double));
                            } else if (extension == ".bgfs") {
                                unsigned int int_val;
                                In.read((char *) (&int_val), sizeof(unsigned int));
                                val = int_val;
                            }
                            graph.ReadNodeFeatures(val, j, graphsNodeFeatureNames[i][k]);

                        }
                    }
                    //Read the edges
                    INDEX added_edges = 0;
                    INDEX original_edges = graph.edges();
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
                        if (graph.ReadEdges(Src, Dst, edgeData)) {
                            ++added_edges;
                        }
                    }
                    graph.set_edges(added_edges);
                    graph.InitLabels();
                }
            }
            In.close();
        }
    }
}

template<typename T>
void GraphData<T>::Load(std::vector<T>& graphs, const std::string &graphPath, GraphFormat graphFormat) {
    switch (graphFormat) {
        case GraphFormat::BGF:
            LoadBGF(graphs, graphPath, graphFormat);
            break;
        case GraphFormat::BGFS:
            LoadBGF(graphs, graphPath, graphFormat);
            break;
        case GraphFormat::BINARY:
            break;
        case GraphFormat::EDGES:
            break;
        case GraphFormat::PEREGRINE_DATA:
            break;
        case GraphFormat::PEREGRINE_SMALL:
            break;
        case GraphFormat::DIMACS:
            break;
        case GraphFormat::AIDS:
            LoadAIDS(graphs, graphPath);
            break;
    }
}

/// Construct new undirected graph
/// \param size
/// \param labels
inline GraphStruct::GraphStruct(INDEX size, const Labels& labels) {
    CreateGraph(size, labels);
}

/// Initialize a new undirected graph
/// \param size
/// \param labels
inline void GraphStruct::CreateGraph(INDEX size, const Labels& labels) {
    this->_nodes = size;
    this->_edges = 0;
    this->_labels = labels;
    this->_graph = std::vector<std::vector<INDEX>>(nodes(), std::vector<INDEX>());
    this->_degrees = std::vector<INDEX>(nodes(), 0);
    this->isTree = false;
}

/// Initialize labeled graph using the label type
/// \param label
inline void GraphStruct::UpdateGraphLabels(LABEL_TYPE label) {
    if (label == LABEL_TYPE::LABELED_SPARSE) {
        this->labeledVectorGraph = std::vector<std::vector<std::vector<INDEX>>>(nodes(), std::vector<std::vector<INDEX>>(_numLabels, std::vector<INDEX>()));
        this->labeledDegreeVector = std::vector<std::vector<INDEX>>(nodes(), std::vector<INDEX>(_numLabels, 0));
    } else if (label == LABEL_TYPE::LABELED_DENSE) {
        this->labeledMapGraph = std::vector<std::unordered_map<Label, Nodes>>(nodes(),std::unordered_map<Label, Nodes>{});
        this->labeledDegreeMap = std::vector<std::unordered_map<Label, INDEX>>(nodes(), std::unordered_map<Label, INDEX>{});
    }
    //set labels of neighbor Nodes
    for (NodeId Id = 0; Id < nodes(); ++Id) {
        for (NodeId neighborId: _graph[Id]) {
            //Create labeled graph
            if (label == LABEL_TYPE::LABELED_DENSE) {
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
            } else if (label == LABEL_TYPE::LABELED_SPARSE) {
                Label NodeLabel = _labels[neighborId];
                this->labeledVectorGraph[Id][NodeLabel].emplace_back(neighborId);
                ++this->labeledDegreeVector[Id][NodeLabel];
            }
        }
    }
}

/// Get the degree of a graph node
/// \param nodeId
/// \return
inline INDEX GraphStruct::degree(NodeId nodeId) const {
    return this->_degrees[nodeId];
}

/// Get the number of graph nodes
/// \return
inline INDEX GraphStruct::nodes() const {
    return _nodes;
}

/// Get the number of graph edges
/// \return
inline INDEX GraphStruct::edges() const {
    return _edges;
}



/// Check if a graph contains some edge
/// \param source
/// \param destination
/// \return
inline bool GraphStruct::edge(NodeId source, NodeId destination) const {
    if (this->degree(source) < this->degree(destination)){
        return std::binary_search(_graph[source].begin(), _graph[source].end(), destination);
    }
    else{
        return std::binary_search(_graph[destination].begin(), _graph[destination].end(), source);
    }
}

/// Get all neighbors of a graph node (const)
/// \param nodeId
/// \return
inline const Nodes &GraphStruct::get_neighbors(NodeId nodeId) const {
    return _graph[nodeId];
}

/// Get all neighbors of a graph node (non const)
/// \param nodeId
/// \return
inline Nodes &GraphStruct::neighbors(NodeId nodeId) {
    return _graph[nodeId];
}

/// Add an edge to an undirected graph
/// \param source
/// \param destination
/// \return
inline bool GraphStruct::add_edge(NodeId source, NodeId destination) {
    if (!edge(source, destination)){
        INDEX ElementId;

        this->_graph[source].emplace_back(destination);
        ElementId = (INDEX) this->_graph[source].size() - 1;
        while (ElementId > 0 && this->_graph[source][ElementId] < this->_graph[source][ElementId - 1]){
            std::swap(this->_graph[source][ElementId], this->_graph[source][ElementId - 1]);
            --ElementId;
        }

        this->_graph[destination].emplace_back(source);
        ElementId = (INDEX) this->_graph[destination].size() - 1;
        while (ElementId > 0 && this->_graph[destination][ElementId] < this->_graph[destination][ElementId - 1]){
            std::swap(this->_graph[destination][ElementId], this->_graph[destination][ElementId - 1]);
            --ElementId;
        }

        ++this->_edges;
        ++this->_degrees[source];
        ++this->_degrees[destination];
        this->maxDegree = std::max(this->maxDegree, (int) std::max(this->degree(source), this->degree(destination)));
        return true;
    }
    return false;
}

void GraphStruct::set_graph(const std::vector<Nodes> &nodes) {
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

/// Get labels of the graph nodes
/// \return
inline const std::vector<Label>& GraphStruct::labels() const {
    return _labels;
}

/// Get label of a graph node
/// \param nodeId
/// \return
inline Label GraphStruct::label(NodeId nodeId) const{
    return _labels[nodeId];
}


/// Get specific neighbor of a graph node
/// \param nodeId
/// \param neighborIdx
/// \return
inline NodeId GraphStruct::neighbor(NodeId nodeId, INDEX neighborIdx) const {
    return _graph[nodeId][neighborIdx];
}

/// Get specific neighbor of a graph node with specific label
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

/// Compare the labeled degrees of graph node with given labeled degrees
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


/// Check if graph node has neighbor with specific label
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

/// Construct graph from path (optionally with node labels)
/// \param graphPath
/// \param relabeling
/// \param withLabels
/// \param labelPath
inline GraphStruct::GraphStruct(const std::string &graphPath,bool relabeling, bool withLabels, const std::string & labelPath, const std::string& format, const std::string& search_name) {
Load(graphPath, relabeling, withLabels, labelPath, format, search_name);
}

/// Save graph in certain path and format
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

bool GraphStruct::ReadEdges(INDEX Src, INDEX Dst, std::vector<double>& edgeData) {
    return this->add_edge(Src, Dst);

}

void GraphStruct::ReadNodeFeatures(double value, int pos, const std::string &nodeFeatureName) {
    if (nodeFeatureName == "label") {
        this->_labels.emplace_back(value);
    }
}

void GraphStruct::Init(const std::string &name, int size, int edges, int nodeFeatures, int edgeFeatures,
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

/// Add labels to a graph
/// \param labels
inline void GraphStruct::set_labels(const Labels* labels) {
    this->_labels = *labels;
    if (this->_labels.size() == this->_graph.size()){
        labelMap = GraphFunctions::GetGraphLabelMap(_labels);
        labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
        this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size()) : this->_numLabels;
        if (this->_numLabels >= 10){
            labelType = LABEL_TYPE::LABELED_DENSE;
        }
        else{
            labelType = LABEL_TYPE::LABELED_SPARSE;
        }
        UpdateGraphLabels(labelType);
    }
    else{
        //TODO throw exception
    }
}

/// Add nodes with labels to a graph
/// \param number
/// \param labels
/// \return
inline NodeId GraphStruct::add_node(INDEX number, Labels* labels) {
    for(INDEX i = 0; i < number; i++) {
        this->_graph.emplace_back();
        this->_degrees.emplace_back(0);
        ++this->_nodes;
        if (labels != nullptr && labels->size() == number){
            this->_labels.emplace_back((*labels)[i]);
        }
    }
    if (labels != nullptr && labels->size() == number){
        this->UpdateGraphLabels(LABEL_TYPE::LABELED_SPARSE);
    }
    return (INDEX) this->_graph.size() - 1;
}

inline bool GraphStruct::CheckSpanningTree(const GraphStruct &spanningTree) const {
    if (spanningTree.nodes() != this->nodes() || !spanningTree.IsTree()) {
        return false;
    }
    NodeId Counter = 0;
    for (auto const & Nodes : spanningTree.graph()) {
        for (auto Node : Nodes) {
            if (!this->edge(Counter, Node)) {
                return false;
            }
        }
        ++Counter;
    }
    return true;
}

inline bool GraphStruct::IsTree() const {
    for (INDEX i = 0; i < this->nodes(); ++i) {
        if(this->degree(i) == 0){
            return false;
        }
    }
    return this->nodes() == this->edges() + 1;
}

inline void GraphStruct::BFSDistances(const GraphStruct &graph, INDEX root, std::vector<INDEX> &distances) {
    distances.clear();
    distances.resize(graph.nodes(), -1);
    std::vector<bool> visitedNodes = std::vector<bool>(graph.nodes(), false);
    visitedNodes[root] = true;
    distances[root] = 0;
    std::deque<NodeId> nodes;
    nodes.push_back(root);
    while (!nodes.empty()){
        NodeId currentNode = nodes.back();
        nodes.pop_back();
        for (NodeId i = 0; i < graph.degree(currentNode); ++i) {
            NodeId neighbor = graph.neighbor(currentNode, i);
            if (!visitedNodes[neighbor]){
                visitedNodes[neighbor] = true;
                distances[neighbor] = distances[currentNode] + 1;
                nodes.push_front(neighbor);
            }
        }
    }
}

inline bool GraphStruct::ReachableNodes(const GraphStruct &graph, NodeId root, std::vector<INDEX> & reachability, INDEX Id, INDEX& number){
    number = 1;
    reachability[root] = Id;
    std::deque<NodeId> nodes;
    nodes.push_back(root);
    while (!nodes.empty()){
        NodeId currentNode = nodes.back();
        nodes.pop_back();
        for (NodeId i = 0; i < graph.degree(currentNode); ++i) {
            NodeId neighbor = graph.neighbor(currentNode, i);
            if (reachability[neighbor] != Id){
                reachability[neighbor] = Id;
                nodes.push_front(neighbor);
                ++number;
            }
        }
    }
    return number == graph.nodes();
}

inline bool GraphStruct::IsConnected(const GraphStruct& graph){
    std::vector<INDEX> reachability = std::vector<INDEX>(graph.nodes(), -1);
    INDEX number;
    return GraphStruct::ReachableNodes(graph, 0, reachability, 0, number);
}

inline void GraphStruct::GetComponents(const GraphStruct& graph, std::vector<INDEX> & components){
    components.clear();
    components.resize(graph.nodes(), -1);
    INDEX number = 0;
    INDEX pos = 0;
    INDEX Id = 0;
    while (number < graph.nodes()){
        INDEX compSize = 0;
        while(components[pos] != -1){
            ++pos;
        }
        GraphStruct::ReachableNodes(graph, pos, components, Id, compSize);
        number += compSize;
        ++Id;
        ++pos;
    }

}

inline void GraphStruct::GetComponents(const GraphStruct& graph, std::vector<Nodes> & components){
    std::vector<INDEX> com;
    GraphStruct::GetComponents(graph, com);
    for (INDEX i = 0; i < com.size(); ++i) {
        INDEX componentId = com[i];
        if (components.empty() || components.size() <= componentId){
            components.resize(componentId + 1, Nodes());
        }
        components[componentId].emplace_back(i);
    }
}

inline void GraphStruct::GetLargestComponent(const GraphStruct& graph, Nodes &nodes) {
    std::vector<Nodes> components;
    GraphStruct::GetComponents(graph, components);
    INDEX maxSize = 0;
    INDEX maxIndex = 0;
    for (INDEX i = 0; i < components.size(); ++i) {
        auto compSize = (INDEX) components[i].size();
        if (compSize > maxSize){
            maxSize = compSize;
            maxIndex = i;
        }
    }
    nodes = components[maxIndex];
}

inline GraphStruct GraphStruct::SubGraph(const GraphStruct& graph, const Nodes& nodeIds){
    GraphStruct g = GraphStruct((INDEX) nodeIds.size(), {});
    std::unordered_map<INDEX, INDEX> idMap = std::unordered_map<INDEX, INDEX>();
    for (INDEX i = 0; i < nodeIds.size(); ++i) {
        idMap[nodeIds[i]] = i;
    }
    for (INDEX srcNode : nodeIds) {
        for (INDEX dstNode : graph.graph()[srcNode]) {
            if (std::find(nodeIds.begin(), nodeIds.end(), dstNode) != nodeIds.end()){
                g.add_edge(idMap[srcNode], idMap[dstNode]);
            }
        }
    }
    return g;
}

inline GraphStruct GraphStruct::GetLargestComponent(const GraphStruct& graph){
    Nodes nodes;
    GraphStruct::GetLargestComponent(graph, nodes);
    return SubGraph(graph, nodes);
}


/// Construct a graph database from a graph path with labels
/// \param graphPath
/// \param labelPath
/// \param extension
/// \param graphSizes
/// \param patternNum
template<typename T>
inline GraphData<T>::GraphData(const std::string& graphPath, const std::string& labelPath, const std::string& searchName, const std::string& extension, bool sort, std::set<int>* graphSizes, int patternNum) {
    std::string ext = extension;
    std::vector<std::string> files;
    if (!searchName.empty()){
        for (const auto &entry: std::filesystem::directory_iterator(graphPath)) {
            if (entry.is_regular_file()) {
                files.emplace_back(entry.path().string());
            }
        }
        for (const auto& path : files) {
            std::string lower_name = searchName;
            std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
                           [](unsigned char c){ return std::tolower(c); } // correct
            );
            std::string path_lower = path;
            std::transform(path_lower.begin(), path_lower.end(), path_lower.begin(),
                           [](unsigned char c){ return std::tolower(c); } // correct
            );
            if (path_lower.find(lower_name) != std::string::npos){
                if (ext.empty()){
                    ext = std::filesystem::path(graphPath).extension().string();
                }
                if (ext == ".bgf" || ext == ".bgfs"){
                    Load(path);
                }else {
                    add(graphPath, labelPath, searchName, ext, sort, graphSizes, patternNum);
                }
            }
        }
    }
    else {
        if (ext.empty()) {
            ext = std::filesystem::path(graphPath).extension().string();
        }
        if (ext == ".bgf" || ext == ".bgfs") {
            Load(graphPath);
        } else {
            add(graphPath, labelPath, searchName, ext, sort, graphSizes, patternNum);
        }
    }
}

/// Construct a graph database from a graph path with labels
/// \param graphPath
/// \param labelPath
/// \param extension
/// \param graphSizes
/// \param patternNum
template<typename T>
inline GraphData<T>::GraphData(GraphFormat graphFormat, const std::string& graphPath, const std::string& search_name, bool sort){
    std::vector<std::string> possible_files;
    std::vector<std::string> files;
    if (!is_regular_file(std::filesystem::path(graphPath))){
        for (const auto &entry: std::filesystem::directory_iterator(graphPath)) {
            if (entry.is_regular_file()) {
                std::string string_path = std::filesystem::path(entry).string();

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
    }
    else{
        files = {graphPath};
    }



    for (auto const & file : files) {
        try {
            Load(this->graphData, file, graphFormat);
        }
        catch(const std::domain_error& e) {
            std::cerr << e.what();
        }
    }
    if (sort) {
        std::sort(this->graphData.begin(), this->graphData.end(), GraphStruct::sort_by_size());
    }
}

/// Add a graph to a graph database
/// \param graph
template<typename T>
inline void GraphData<T>::add(const T& graph) {
    graphData.push_back(graph);
}

/// Add graphs from path to a graph database, optionally with labels
/// \param graphPath
/// \param labelPath
/// \param extension
/// \param graphSizes
/// \param patternNum
template<typename T>
inline void
GraphData<T>::add(const std::string &graphPath, const std::string &labelPath, const std::string& searchName, const std::string & extension, bool sort, std::set<int> *graphSizes, int patternNum) {
    GraphStruct::LoadGraphsFromPath(graphPath, labelPath, this->graphData, searchName, extension, sort, graphSizes, patternNum);
}

/// Add a vector of undirected graph to a graph
/// \param graphs
template<typename T>
inline void GraphData<T>::add(const std::vector<T> &graphs) {
    this->graphData.insert(this->graphData.end(), graphs.begin(), graphs.end());
}

/// Add graph from path a graph database
/// \param graphPath
/// \param withLabels
template<typename T>
inline void GraphData<T>::add_graph(const std::string& graphPath, bool withLabels) {
    this->graphData.emplace_back(graphPath, withLabels);
}

template<typename T>
inline void GraphStruct::LoadGraphsFromPath(const std::string& graphPath, const std::string& labelPath, std::vector<T> &graphs, const std::string & searchName, const std::string & extension, const bool sort, std::set<int>* graphSizes, int patternNum)
{
    std::unordered_map<INDEX, INDEX> sizesNumMap;
    if (graphSizes != nullptr) {
        for (INDEX size : *graphSizes) {
            sizesNumMap.insert({size, 0});
        }
    }
        for (const auto &entry: std::filesystem::directory_iterator(graphPath)) {
            if (std::filesystem::is_regular_file(entry)) {
                if (entry.path().extension().string() == extension || extension.empty()) {
                    std::string path = entry.path().string();
                    if (searchName.empty() || path.find(searchName) != std::string::npos) {
                        T graph(path, true, true, labelPath);
                        if (graphSizes == nullptr ||
                            graphSizes->find(static_cast<int>(graph.nodes())) != graphSizes->end()) {
                            if (patternNum == -1 || graphSizes == nullptr ||
                                sizesNumMap[graph.nodes()] < patternNum) {
                                graphs.emplace_back(graph);
                                if (graphSizes != nullptr) {
                                    ++sizesNumMap[graph.nodes()];
                                }
                            }
                        }
                    }
                }
            }
        }
        if (sort) {
            std::sort(graphs.begin(), graphs.end(), GraphStruct::sort_by_name());
        }

}

void GraphStruct::Convert(const std::string &graphPath, const std::string &Name,
                          GraphFormat Format, bool relabeling, bool withLabels, const std::string &labelPath, bool Labeled) {
    GraphStruct graphStruct = GraphStruct(graphPath, relabeling, withLabels, labelPath);
    graphStruct.Save({"", Name, Format, Labeled});
}

void GraphStruct::Convert(const std::string &path, GraphFormat Format, const std::string &extension) {
    std::vector<std::string> files;
    for (const auto &entry: std::filesystem::directory_iterator(path)) {
        if (entry.is_regular_file()) {
            files.emplace_back(entry.path().string());
        }
    }
    for(const auto& f_path : files){
        if (std::filesystem::path(f_path).extension() == extension){
            GraphStruct graphStruct = GraphStruct(f_path, true, false, "");
            graphStruct.Save({"", "", Format, false});
        }
    }

}

void GraphStruct::Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format, const std::string& search_name) {
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

            if (GraphStruct::ReadBGF(extension, In, saveVersion, graphNumber, graphsNames, graphsTypes, graphsSizes,
                                     graphsNodeFeatureNames, graphsEdges, graphsEdgeFeatureNames)) {
                for (int i = 0; i < graphNumber; ++i) {
                    //Create graph
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
                            if (this->ReadEdges(Src, Dst, edgeData)) {
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
                    labelMap = GraphFunctions::GetGraphLabelMap(_labels);
                    labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
                    this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size()) : this->_numLabels;
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
                if (!GraphStruct::add_edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (this->edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base graph)!"
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
                    labelMap = GraphFunctions::GetGraphLabelMap(_labels);
                    labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
                    this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size()) : this->_numLabels;
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
                if (!GraphStruct::add_edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (this->edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            if (num_edges_duplicates > 0) {
                std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base graph)!"
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
                    labelMap = GraphFunctions::GetGraphLabelMap(_labels);
                    labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
                    this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size())
                                                                : this->_numLabels;
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
                if (!GraphStruct::add_edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (this->edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            if (num_edges_duplicates > 0) {
                std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base graph)!"
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
                    labelMap = GraphFunctions::GetGraphLabelMap(_labels);
                    labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
                    this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size())
                                                                : this->_numLabels;
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
    }

    this->sortNeighborIds();
    this->isTree = IsTree();

}

/// Construct a directed graph from an undirected graph
/// \param graph
inline DGraphStruct::DGraphStruct(GraphStruct & graph) : GraphStruct(graph){
    for (INDEX i = 0; i < graph.nodes(); ++i) {
        this->_in_degrees.emplace_back(graph.degree(i));
        this->_out_degrees.emplace_back(graph.degree(i));
    }
}

/// Construct a directed graph by certain size
/// \param size
inline DGraphStruct::DGraphStruct(INDEX size) : GraphStruct(size, Labels()){
    this->_in_degrees.resize(size, 0);
    this->_out_degrees.resize(size, 0);
}

/// Read directed graph from file
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

/// Constructor of directed graph with certain size and given labels
/// \param size
/// \param labels
inline DGraphStruct::DGraphStruct(INDEX size, const Labels &labels) : GraphStruct(size, labels) {
    this->_in_degrees.resize(size, 0);
    this->_out_degrees.resize(size, 0);
}

/// Ad an edge to a directed graph
/// \param source
/// \param destination
/// \return
inline bool DGraphStruct::add_edge(NodeId source, NodeId destination) {
    if (!edge(source, destination)){
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

/// Check if edge in directed graph exists
/// \param source
/// \param destination
/// \return
inline bool DGraphStruct::edge(NodeId source, NodeId destination) const {
    return DGraphStruct::edge(source, destination, true);
}

/// Check if edge in directed graph exists
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
/// \param rootNodeId
/// \return
inline DGraphStruct DGraphStruct::GetBFSTree(const GraphStruct &graph, NodeId rootNodeId) {
        DGraphStruct BFSTree = DGraphStruct(graph.nodes(), graph.labels());
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
                    BFSTree.add_edge(NextNodeId, neighborId);
                }
            }
        }
        return BFSTree;
}

void DGraphStruct::Load(const std::string &graphPath, bool relabeling, bool withLabels, const std::string &labelPath, const std::string& format) {
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
                    labelMap = GraphFunctions::GetGraphLabelMap(_labels);
                    labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
                    this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size()) : this->_numLabels;
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
                GraphStruct::add_edge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second]);
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
                    labelMap = GraphFunctions::GetGraphLabelMap(_labels);
                    labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
                    this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size()) : this->_numLabels;
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
    }

    this->sortNeighborIds();
    this->isTree = IsTree();
}

void DGraphStruct::Save(const SaveParams& saveParams) {
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

void DGraphStruct::Init(const std::string &name, int size, int edges, int nodeFeatures, int edgeFeatures,
                        const std::vector<std::string> &nodeFeatureNames,
                        const std::vector<std::string> &edgeFeatureNames) {
    GraphStruct::Init(name, size, edges, nodeFeatures, edgeFeatures, nodeFeatureNames, edgeFeatureNames);
    //Create graph
    this->_in_degrees.resize(_nodes);
    this->_out_degrees.resize(_nodes);
}

bool DGraphStruct::ReadEdges(INDEX Src, INDEX Dst, std::vector<double> &edgeData) {
    return this->add_edge(Src, Dst);
}

void DGraphStruct::ReadNodeFeatures(double value, int pos, const std::string &nodeFeatureName) {
    GraphStruct::ReadNodeFeatures(value, pos, nodeFeatureName);
}

void DGraphStruct::WriteEdges(std::ofstream& Out, const SaveParams& saveParams){
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

/// Default constructor
inline DDataGraph::DDataGraph() = default;

/// Read directed data graph from file TODO binary extension
/// \param graphPath path of the graph file
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
        this->isTree = IsTree();
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

/// Get edge data for some edge in a directed data graph
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

/// Add an data edge to a directed data graph
/// \param source
/// \param destination
/// \param data
/// \return
inline bool DDataGraph::add_edge(NodeId source, NodeId destination, std::vector<double>& data) {
    add_edge_data(EDGE{source, destination}, data);
    return DGraphStruct::add_edge(source, destination);
}


/// Get the node data for some node in the directed data graph (by type)
/// \param node
/// \param type
/// \return
inline double DDataGraph::get_node_data(NodeId node, const std::string &type) const {
    return get_node_data(node, _node_data_names.at(type));
}

/// Get the node data for some node in the directed data graph (by index)
/// \param node
/// \param index
/// \return
inline double DDataGraph::get_node_data(NodeId node, int index) const {
    return _node_data[node][index];
}

/// Get the node data for some node in the directed data graph (all data)
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

void DDataGraph::Init(const std::string &name, int size, int edges, int nodeFeatures, int edgeFeatures,
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

void DDataGraph::ReadNodeFeatures(double value, int pos, const std::string &nodeFeatureName) {
    this->_node_data[pos].emplace_back(value);
    if (nodeFeatureName == "label") {
        this->_labels.emplace_back((int) value);
    }
}

bool DDataGraph::ReadEdges(INDEX Src, INDEX Dst, std::vector<double> &edgeData) {
    return this->add_edge(Src, Dst, edgeData);
}


void DDataGraph::WriteGraph(std::ofstream& Out, const SaveParams& saveParams){
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

    Out.write((char *) (&numFeatures), sizeof(unsigned int));
    for (int j = 0; j < numFeatures; ++j) {
        unsigned int edgeFeatureStringLength = featureNames[j].length();
        Out.write((char *) (&edgeFeatureStringLength), sizeof(edgeFeatureStringLength));
        Out.write(featureNames[j].c_str(), edgeFeatureStringLength);
    }
}

void DDataGraph::WriteNodeFeatures(std::ofstream& Out,const SaveParams& saveParams){
    for (int j = 0; j < this->nodes(); ++j) {
        for (int k = 0; k < _node_data_size; ++k) {
            auto val = _node_data[j][k];
            Out.write((char *) (&val), sizeof(double));
        }
    }
}
void DDataGraph::WriteEdges(std::ofstream& Out,const SaveParams& saveParams){
    INDEX Src = 0;
    for (auto const &edges: this->_graph) {
        for (auto Dst: edges) {
            Out.write((char *) (&Src), sizeof(INDEX));
            Out.write((char *) (&Dst), sizeof(INDEX));
            for (int j = 0; j < _edge_data_size; ++j) {
                auto val = _edge_data[Src][Dst][j];
                Out.write((char *) (&val), sizeof(INDEX));
            }
        }
        ++Src;
    }
}

/// Deprecated Save graph in certain path and format
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

bool GraphStruct::ReadBGF(const std::string& extension,std::ifstream& In, int &saveVersion,  int &graphNumber, std::vector<std::string> &graphsNames,
                          std::vector<GraphType> &graphsTypes, std::vector<INDEX> &graphsSizes,
                          std::vector<std::vector<std::string>> &graphsNodeFeatureNames,
                          std::vector<INDEX> &graphsEdges,
                          std::vector<std::vector<std::string>> &graphsEdgeFeatureNames) {

    In.read((char *) (&saveVersion), sizeof(int));
    if (saveVersion == 1) {
        In.read((char *) (&graphNumber), sizeof(int));
        for (int i = 0; i < graphNumber; ++i) {
            std::string graphName;
            unsigned int stringLength;
            GraphType Type;
            INDEX Size;
            unsigned int nodeFeatures;
            std::vector<std::string> nodeFeatureNames;

            In.read((char *) (&stringLength), sizeof(unsigned int));
            graphName.resize(stringLength);
            In.read((char *) graphName.c_str(), stringLength);
            graphsNames.emplace_back(graphName);
            In.read((char *) (&Type), sizeof(GraphType));
            graphsTypes.emplace_back(Type);
            if (extension == ".bgf") {
                In.read((char *) (&Size), sizeof(INDEX));
            }
            else if (extension == ".bgfs"){
                unsigned int int_size;
                In.read((char *) (&int_size), sizeof(unsigned int));
                Size = int_size;
            }
            graphsSizes.emplace_back(Size);

            In.read((char *) (&nodeFeatures), sizeof(unsigned int));
            for (int j = 0; j < nodeFeatures; ++j) {
                unsigned int nodeFeatureStringLength;
                std::string nodeFeatureName;
                In.read((char *) (&nodeFeatureStringLength), sizeof(unsigned int));
                nodeFeatureName.resize(nodeFeatureStringLength);
                In.read((char *) nodeFeatureName.c_str(), nodeFeatureStringLength);
                nodeFeatureNames.emplace_back(nodeFeatureName);
            }
            graphsNodeFeatureNames.emplace_back(nodeFeatureNames);

            INDEX Edges;
            if (extension == ".bgf") {
                In.read((char *) (&Edges), sizeof(INDEX));
            }
            else if (extension == ".bgfs"){
                unsigned int int_edges;
                In.read((char *) (&int_edges), sizeof(unsigned int));
                Edges = int_edges;
            }
            graphsEdges.emplace_back(Edges);
            unsigned int edgeFeatures;
            std::vector<std::string> edgeFeatureNames;
            In.read((char *) (&edgeFeatures), sizeof(unsigned int));
            for (int j = 0; j < edgeFeatures; ++j) {
                unsigned int edgeFeatureStringLength;
                std::string edgeFeatureName;
                In.read((char *) (&edgeFeatureStringLength), sizeof(unsigned int));
                edgeFeatureName.resize(edgeFeatureStringLength);
                In.read((char *) edgeFeatureName.c_str(), edgeFeatureStringLength);
                edgeFeatureNames.emplace_back(edgeFeatureName);
            }
            graphsEdgeFeatureNames.emplace_back(edgeFeatureNames);
        }
        return true;
    }
    return false;
}

void GraphStruct::InitLabels() {
    //Add graph labels if possible
    if (this->_labels.size() == this->_graph.size()) {
        labelMap = GraphFunctions::GetGraphLabelMap(_labels);
        labelFrequencyMap = GraphFunctions::GetLabelFrequency(labelMap);
        this->_numLabels = (this->_numLabels == -1) ? static_cast<int>(labelMap.size()) : this->_numLabels;
        if (this->_numLabels >= 10) {
            labelType = LABEL_TYPE::LABELED_DENSE;
        } else {
            labelType = LABEL_TYPE::LABELED_SPARSE;
        }
        UpdateGraphLabels(labelType);
    } else {
        this->labelType = LABEL_TYPE::UNLABELED;
    }
}

void GraphStruct::WriteGraph(std::ofstream& Out, const SaveParams& saveParams) {
    const std::string graphName = this->_name;
    unsigned int stringLength = graphName.length();
    GraphType Type = this->graphType;
    auto Size = this->nodes();
    unsigned int nodeFeatures;
    std::vector<std::string> nodeFeatureNames;

    Out.write((char *) (&stringLength), sizeof(unsigned int));
    Out.write(graphName.c_str(), stringLength);
    Out.write((char *) (&Type), sizeof(GraphType));

    if(saveParams.Format == GraphFormat::BGF) {
        Out.write((char *) (&Size), sizeof(INDEX));
    }
    else if(saveParams.Format == GraphFormat::BGFS) {
        auto int_size = (unsigned int) Size;
        Out.write((char *) (&int_size), sizeof(unsigned int));
    }

    nodeFeatures = 0;
    if (saveParams.Labeled || _labels.size() == _nodes) {
        nodeFeatures = 1;
        nodeFeatureNames.emplace_back("label");
    }
    Out.write((char *) (&nodeFeatures), sizeof(unsigned int));

    for (int j = 0; j < nodeFeatures; ++j) {
        unsigned int nodeFeatureStringLength = nodeFeatureNames[j].length();
        Out.write((char *) (&nodeFeatureStringLength), sizeof(unsigned int));
        Out.write(nodeFeatureNames[j].c_str(), nodeFeatureStringLength);
    }

    auto Edges = this->edges();
    if(saveParams.Format == GraphFormat::BGF) {
        Out.write((char *) (&Edges), sizeof(INDEX));
    }
    else if(saveParams.Format == GraphFormat::BGFS) {
        auto int_edges = (unsigned int) Edges;
        Out.write((char *) (&int_edges), sizeof(unsigned int));
    }

    unsigned int edgeFeatures = 0;
    std::vector<std::string> edgeFeatureNames = {};
    Out.write((char *) (&edgeFeatures), sizeof(unsigned int));
    for (int j = 0; j < edgeFeatures; ++j) {
        unsigned int edgeFeatureStringLength = edgeFeatureNames[j].length();
        Out.write((char *) (&edgeFeatureStringLength), sizeof(unsigned int));
        Out.write(edgeFeatureNames[j].c_str(), edgeFeatureStringLength);
    }

}

void GraphStruct::WriteNodeFeatures(std::ofstream& Out, const SaveParams &saveParams) {
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

void GraphStruct::WriteEdges(std::ofstream& Out, const SaveParams &saveParams) {
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

void GraphStruct::DFS(const GraphStruct &graph, GraphStruct &tree, Nodes &nodeOrder, NodeId rootNodeId, int seed) {
    std::mt19937_64 gen(seed);
    NodeId max = -1;
    if (rootNodeId == max) {
        rootNodeId = std::uniform_int_distribution<INDEX>(0, graph.nodes() - 1)(gen);
    }
    std::deque<std::pair<INDEX, INDEX>> swappedIds;
    std::vector<NodeId> neighborIds = std::vector<NodeId>(graph.maxDegree, 0);
    std::iota(neighborIds.begin(), neighborIds.end(), 0);
    for (int i = 0; i < graph.nodes(); ++i) {
        tree.add_node();
    }
    tree.graphType = GraphType::TREE;
    std::unordered_set<INDEX> visitedNodes;
    std::vector<NodeId> CurrentNodes;
    CurrentNodes.emplace_back(rootNodeId);
    nodeOrder.emplace_back(rootNodeId);
    while (!CurrentNodes.empty()) {
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        if (visitedNodes.find(NextNodeId) == visitedNodes.end()) {
            visitedNodes.insert(NextNodeId);
            nodeOrder.emplace_back(NextNodeId);
            INDEX degree = graph.get_neighbors(NextNodeId).size();
            for (INDEX i = 0; i < degree; ++i) {
                //Get random neighbor
                INDEX idx = std::uniform_int_distribution<INDEX>(i, degree - 1)(gen);
                NodeId NeighborNodeId = graph.get_neighbors(NextNodeId)[neighborIds[idx]];
                std::swap(neighborIds[idx], neighborIds[i]);
                swappedIds.push_front(std::pair<int, int>{idx, i});
                if (visitedNodes.find(NeighborNodeId) == visitedNodes.end()) {
                    CurrentNodes.emplace_back(NeighborNodeId);
                }
            }
            tree.add_edge(NextNodeId, CurrentNodes.back());
            for (auto const &[a, b]: swappedIds) {
                std::swap(neighborIds[b], neighborIds[a]);
            }
            swappedIds.clear();
        }
    }
}

void GraphStruct::OptOrdering(const GraphStruct &graph, Nodes &nodeOrder) {
    NodeId rootNodeId;
    NodeId max_degree = 0;
    for (int i = 0; i < graph.nodes(); ++i) {
        if (graph.degree(i) > max_degree){
            max_degree = graph.degree(i);
            rootNodeId = i;
        }
    }
    std::unordered_set<INDEX> visitedNodes = {rootNodeId};
    std::vector<NodeId> CurrentNodes;
    CurrentNodes.emplace_back(rootNodeId);
    nodeOrder.emplace_back(rootNodeId);
    while (!CurrentNodes.empty()){
        NodeId NextNodeId = CurrentNodes.back();
        CurrentNodes.pop_back();
        std::vector<std::pair<NodeId, INDEX>> node_degrees;
        for(auto n : graph.get_neighbors(NextNodeId)){
            node_degrees.emplace_back(n, graph.degree(n));
        }
        std::sort(node_degrees.begin(), node_degrees.end(), [&](const auto& a, const auto& b)
        {
            return a.second > b.second;
        });
        for (auto [n, d] : node_degrees) {
            if (visitedNodes.find(n) == visitedNodes.end()){
                CurrentNodes.emplace_back(n);
                nodeOrder.emplace_back(n);
                visitedNodes.insert(n);
            }
        }
    }
}

void GraphStruct::ReorderGraph(GraphStruct &graph, const Nodes &nodeOrder) {

    std::unordered_map<NodeId, NodeId> nodeMap;
    int counter = 0;
    for (auto node : nodeOrder) {
        nodeMap.insert({node, counter});
        ++counter;
    }
    std::vector<Nodes> newNodes = std::vector<Nodes>(graph.nodes());
    counter = 0;
    for (auto node : nodeOrder) {
        const Nodes& neighbors = graph.get_neighbors(node);
        for (auto n : neighbors) {
            newNodes[counter].emplace_back(nodeMap[n]);
        }
        ++counter;
    }
    graph.set_graph(newNodes);

}




#endif //HOPS_DATACLASSES_H
