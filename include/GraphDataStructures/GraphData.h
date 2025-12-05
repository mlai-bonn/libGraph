//
// Created by florian on 13.09.25.
//

#ifndef GRAPH_DATA_H
#define GRAPH_DATA_H

#include "GraphStructs.h"

enum class GraphFormat;

template <typename T>
class GraphData {
public:
    GraphData()= default;
    explicit GraphData(const std::string& graphPath, const std::string& labelPath = "", const std::string& searchName = "", const std::string& extension = "", bool sort = true, std::set<int>* graphSizes = nullptr, int patternNum = -1);
    explicit GraphData(GraphFormat graphFormat, const std::string& graphPath, const std::string& search_name = "", bool sort = true);

    void Load(const std::string &graphPath);
    void Load(std::vector<T>& graphs, const std::string &graphPath, GraphFormat graphFormat = GraphFormat::BGFS, const std::string& datasetName = "");

    void LoadBGF(std::vector<T>& graphs, const std::string &graphPath, GraphFormat graphFormat);
    void LoadAIDS(std::vector<T>& graphs, const std::string &graphPath);
    void LoadGRAPHLIST(std::vector<T>& graphs, const std::string &graphPath, const std::string& datasetName = "");
    void Save(const SaveParams& saveParams = SaveParams());

    void add_graph(const std::string& graphPath, bool withLabels = false);
    void add(const T& graph);
    void add(const std::vector<T>& graphs);
    void add(const std::string& graphPath, const std::string& labelPath = "", const std::string& searchName = "", const std::string& extension = ".edges", bool sort = true, std::set<int>* graphSizes = nullptr, int patternNum = -1);

    std::vector<T> graphData;
    T& operator[](size_t index){return graphData[index];};
    [[nodiscard]] INDEX size() const{return static_cast<INDEX>(graphData.size());};

    // Setter and getter
    void SetName(const std::string& name){_name = name;};
    std::string GetName() const{return _name;};

private:
    std::string _name;
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

            if (LoadSave::ReadBGF(extension, In, saveVersion, graphNumber, graphsNames, graphsTypes, graphsSizes,
                                     graphsNodeFeatureNames, graphsEdges, graphsEdgeFeatureNames)) {
                for (int i = 0; i < graphNumber; ++i) {
                    graphs.emplace_back();
                    T& graph = graphs.back();
                    //Create _graph
                    graph.Init(graphsNames[i], (int) graphsSizes[i], (int) graphsEdges[i],
                                   (int) graphsNodeFeatureNames[i].size(), (int) graphsEdgeFeatureNames[i].size(),
                                   graphsNodeFeatureNames[i], graphsEdgeFeatureNames[i]);
                    graph.SetPath(graphPath);
                    //Read the nodes
                    for (int j = 0; j < graphsSizes[i]; ++j) {
                        //std::cout << "Reading node " << j << " of graph " << i << std::endl;
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
                            //std::cout << "Feature values: " << val;
                        }
                        //std::cout << std::endl;
                    }
                    // print nodes loaded for graph i
                    //std::cout << "Loaded " << graphsSizes[i] << " nodes for graph " << i << std::endl;
                    //std::cout << "Reading edges for graph " << i << std::endl;
                    //Read the edges
                    INDEX added_edges = 0;
                    INDEX original_edges = graph.edges();
                    for (INDEX j = 0; j < original_edges; ++j) {
                        INDEX Src = 0;
                        INDEX Dst = 0;
                        if (extension == ".bgf") {
                            In.read(reinterpret_cast<char *>(&Src), sizeof(INDEX));
                            In.read(reinterpret_cast<char *>(&Dst), sizeof(INDEX));
                            // Check for errors
                        } else if (extension == ".bgfs") {
                            unsigned int int_Src;
                            unsigned int int_Dst;
                            In.read(reinterpret_cast<char *>(&int_Src), sizeof(unsigned int));
                            In.read(reinterpret_cast<char *>(&int_Dst), sizeof(unsigned int));
                            Src = (INDEX) int_Src;
                            Dst = (INDEX) int_Dst;
                        }
                        //std::cout << Src << " " << Dst;
                        //if (graphsEdgeFeatureNames.size() > 0) {
                        //    std::cout << " with features: ";
                        //}
                        std::vector<double> edgeData;
                        for (int k = 0; k < graphsEdgeFeatureNames[i].size(); ++k) {
                            double val;
                            In.read((char *) (&val), sizeof(val));
                            edgeData.emplace_back(val);
                            //std::cout << val << " ";
                        }
                        //std::cout << std::endl;
                        if (graph.AddEdge(Src, Dst, edgeData)) {
                            ++added_edges;
                        }
                    }
                    //std::cout << "Graph " << graph.GetName() << " with id " << i << " with " << added_edges << " edges loaded." << std::endl;
                    graph.set_edge_num(added_edges);
                    graph.InitLabels();
                }
            }
    In.close();

}

template<typename T>
void GraphData<T>::LoadGRAPHLIST(std::vector<T> &graphs, const std::string &graphPath, const std::string& datasetName) {
    // Format: First line: number of graphs
    // For each graph:
    // Line with number of nodes
    // Line with node labels (space separated)
    // Lines with edges: nodeA nodeB edgeLabel (space separated)
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
    int graph_counter = 0;
    std::getline(infile, line);
    int total_graphs = std::stoi(line);
    graphs.reserve(total_graphs);
    while (graphs.size() < total_graphs) {
        // show progress every 10 graphs
        if (graph_counter % 100 == 0) {
            std::cout << "Loading graph " << graph_counter + 1 << " / " << total_graphs << std::endl;
        }
        graphEdges.clear();
        nodeLabels.clear();
        edgeLabels.clear();
        graphNodeIds.clear();
        originalIdsToNodeIds.clear();
        graphs.emplace_back();
        T& graph = graphs.back();
        if (num_nodes == 0) {
            // Read number of nodes
            std::getline(infile, line);
            num_nodes = std::stoi(line);
        }
        // Read node labels
        std::getline(infile, line);
        std::istringstream iss_labels(line);
        for (unsigned int i = 0; i < num_nodes; ++i) {
            iss_labels >> a;
            nodeLabels.emplace_back(std::stoi(a));
        }
        graph.Init(datasetName + "_" + std::to_string(graph_counter), 0, 0, 1, 1, {"label"}, {"label"});
        std::vector<std::vector<double>> nodeData;
        for (unsigned int i = 0; i < num_nodes; ++i) {
            nodeData.emplace_back(std::vector<double>{static_cast<double>(nodeLabels[i])});
        }
        graph.AddNodes(num_nodes, nodeLabels, nodeData);
        graph.SetPath(graphPath);
        // Read edges until the line has only one entry or is empty
        num_edges = 0;
        while (std::getline(infile, line) && !line.empty()) {
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            while (iss >> a) {
                tokens.emplace_back(a);
            }
            if (tokens.size() == 1) {
                // Next graph starts, create current graph
                // add edges
                int edgeCounter = 0;
                unsigned int added_edges = 0;
                for (const auto & edge: graphEdges) {
                    if (graph.AddEdge(edge.first, edge.second, {edgeLabels[edgeCounter]})){
                        ++added_edges;
                    }
                    ++edgeCounter;
                }
                graph.set_edge_num(added_edges);
                graph.InitLabels();
                num_nodes = std::stoi(tokens[0]);
                ++graph_counter;
                break;
            }
            graphEdges.emplace_back(std::stoi(tokens[0]), std::stoi(tokens[1]));
            graphNodeIds.emplace(std::stoi(tokens[0]));
            graphNodeIds.emplace(std::stoi(tokens[1]));
            edgeLabels.emplace_back(std::stoi(tokens[2]));
            ++num_edges;
        }

        // add edges
    }
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
            //Create _graph
            graph.Init(std::to_string(graph_counter), num_nodes, 0,1, 1,{"node_label"}, {"edge_label"});

            INDEX nodeCounter = 0;
            for (auto x: graphNodeIds) {
                originalIdsToNodeIds.insert({x, nodeCounter});
                ++nodeCounter;
            }
            unsigned int num_edges_duplicates = 0;
            for (auto edge: graphEdges) {
                if (!graph.AddEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second], {})) {
                    std::cout << "Edge: " << originalIdsToNodeIds[edge.first] << " "
                              << originalIdsToNodeIds[edge.second] << " has not been added because: " << std::endl;
                    if (graph.IsEdge(originalIdsToNodeIds[edge.first], originalIdsToNodeIds[edge.second])) {
                        std::cout << "It already exists!" << std::endl;
                        ++num_edges_duplicates;
                    }
                }
            }
            if (num_edges_duplicates > 0) {
                std::cout << num_edges_duplicates << " edges are not added because of duplicates (directed base _graph)!"
                          << std::endl;
                std::cout << graphEdges.size() - num_edges_duplicates << " edges are remaining" << std::endl;
            }

            int max_label = *std::max(nodeLabels.begin(), nodeLabels.end());
            graph.SetLabels(&nodeLabels);
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
    // remove the .extension from the path if it exists (do not delete relative path .. !!)
    std::string::size_type first_pos = saveName.find_first_of('.');
    std::string::size_type last_pos = saveName.find_last_of('.');
    if (last_pos != std::string::npos && last_pos-first_pos != 1)
    {
        saveName = saveName.substr(0, last_pos);
    }
    int counter;
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
            // iterate over all graphs and save them as list of edges
            counter = 0;
            for (auto &graph: this->graphData) {
                auto out_name = saveName;
                // add _id_counter
                out_name.append("_id_" + std::to_string(counter) + ".edges");
                std::ofstream Out(out_name, std::ios::out);
                for (auto beginEdge = graph.first_edge(); beginEdge != graph.last_edge(); ++beginEdge) {
                    Out << (*beginEdge).first << " " << (*beginEdge).second << std::endl;
                }
                Out.close();
                ++counter;
            }
            break;
        case GraphFormat::PEREGRINE_DATA:
            break;
        case GraphFormat::PEREGRINE_SMALL:
            break;
        case GraphFormat::DIMACS:
            break;
        case GraphFormat::AIDS:
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

            if (LoadSave::ReadBGF(extension, In, saveVersion, graphNumber, graphsNames, graphsTypes, graphsSizes,
                                     graphsNodeFeatureNames, graphsEdges, graphsEdgeFeatureNames)) {
                for (int i = 0; i < graphNumber; ++i) {
                    this->graphData.emplace_back();
                    T& graph = this->graphData.back();
                    //Create _graph
                    graph.Init(graphsNames[i], (int) graphsSizes[i], (int) graphsEdges[i],
                                   (int) graphsNodeFeatureNames[i].size(), (int) graphsEdgeFeatureNames[i].size(),
                                   graphsNodeFeatureNames[i], graphsEdgeFeatureNames[i]);
                    graph.SetPath(graphPath);
                    //Read the nodes
                    for (int j = 0; j < graphsSizes[i]; ++j) {
                        //std::cout << "Reading node " << j << " of graph " << i << std::endl;
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
                            //std::cout << "Feature values: " << val;
                        }
                        //std::cout << std::endl;
                    }
                    // print nodes loaded for graph i
                    //std::cout << "Loaded " << graphsSizes[i] << " nodes for graph " << i << std::endl;
                    //std::cout << "Reading edges for graph " << i << std::endl;
                    //Read the edges
                    INDEX added_edges = 0;
                    INDEX original_edges = graph.edges();
                    for (INDEX j = 0; j < original_edges; ++j) {
                        INDEX Src = 0;
                        INDEX Dst = 0;
                        if (extension == ".bgf") {
                            In.read(reinterpret_cast<char *>(&Src), sizeof(INDEX));
                            In.read(reinterpret_cast<char *>(&Dst), sizeof(INDEX));
                            // Check for errors
                        } else if (extension == ".bgfs") {
                            unsigned int int_Src;
                            unsigned int int_Dst;
                            In.read(reinterpret_cast<char *>(&int_Src), sizeof(unsigned int));
                            In.read(reinterpret_cast<char *>(&int_Dst), sizeof(unsigned int));
                            Src = (INDEX) int_Src;
                            Dst = (INDEX) int_Dst;
                        }
                        //std::cout << Src << " " << Dst;
                        //if (graphsEdgeFeatureNames.size() > 0) {
                        //    std::cout << " with features: ";
                        //}
                        std::vector<double> edgeData;
                        for (int k = 0; k < graphsEdgeFeatureNames[i].size(); ++k) {
                            double val;
                            In.read((char *) (&val), sizeof(val));
                            edgeData.emplace_back(val);
                            //std::cout << val << " ";
                        }
                        //std::cout << std::endl;
                        if (graph.AddEdge(Src, Dst, edgeData)) {
                            ++added_edges;
                        }
                    }
                    //std::cout << "Graph " << graph.GetName() << " with id " << i << " with " << added_edges << " edges loaded." << std::endl;
                    graph.set_edge_num(added_edges);
                    graph.InitLabels();
                }
            }
            In.close();
        }
    }
}

template<typename T>
void GraphData<T>::Load(std::vector<T>& graphs, const std::string &graphPath, GraphFormat graphFormat, const std::string& datasetName) {
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
        case GraphFormat::GRAPHLIST:
            LoadGRAPHLIST(graphs, graphPath, datasetName);
            break;
        case GraphFormat::AIDS:
            LoadAIDS(graphs, graphPath);
            break;
    }
}

/// Construct a _graph database from a _graph path with labels
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

/// Construct a _graph database from a _graph path with labels
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

/// Add a _graph to a _graph database
/// \param graph
template<typename T>
inline void GraphData<T>::add(const T& graph) {
    graphData.push_back(graph);
}

/// Add graphs from path to a _graph database, optionally with labels
/// \param graphPath
/// \param labelPath
/// \param extension
/// \param graphSizes
/// \param patternNum
template<typename T>
inline void
GraphData<T>::add(const std::string &graphPath, const std::string &labelPath, const std::string& searchName, const std::string & extension, bool sort, std::set<int> *graphSizes, int patternNum) {
    LoadSave::LoadGraphsFromPath(graphPath, labelPath, this->graphData, searchName, extension, sort, graphSizes, patternNum);
}

/// Add a vector of undirected _graph to a _graph
/// \param graphs
template<typename T>
inline void GraphData<T>::add(const std::vector<T> &graphs) {
    this->graphData.insert(this->graphData.end(), graphs.begin(), graphs.end());
}

/// Add _graph from path a _graph database
/// \param graphPath
/// \param withLabels
template<typename T>
inline void GraphData<T>::add_graph(const std::string& graphPath, bool withLabels) {
    this->graphData.emplace_back(graphPath, withLabels);
}


#endif //GRAPH_DATA_H