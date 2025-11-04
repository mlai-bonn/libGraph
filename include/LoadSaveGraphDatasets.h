//
// Created by Florian on 16.03.2021.
//

#ifndef LOADSAVEGRAPHDATASETS_H
#define LOADSAVEGRAPHDATASETS_H

#include <fstream>
#include <string>
#include <set>
#include <unordered_map>
#include "typedefs.h"
#include "GraphDataStructures/GraphBase.h"
#include "GraphDataStructures/GraphData.h"
#include "GraphDataStructures/GraphStructs.h"
#include "GraphDataStructures/GraphUndirectedFeaturedBase.h"


class LoadSaveGraphDatasets {
public:

    template <typename T>
    static void LoadTUDortmundGraphData(const std::string& path,
        const std::string& dbName, GraphData<T>& graphs,
        std::vector<int>& graphLabels,
        std::vector<std::vector<Label>>* graphsNodeLabels=nullptr,
        std::vector<std::vector<Label>>* graphsEdgeLabels= nullptr,
        std::vector<std::vector<double>>* graphsNodeAttributes= nullptr,
        std::vector<std::vector<double>>* graphsEdgeAttributes= nullptr);
    static bool PreprocessTUDortmundGraphData(const std::string& dataset_name, const std::string &input_path, const std::string &output_path);
    template <typename T>
    static void LoadPreprocessedTUDortmundGraphData(const std::string& dataset_name , const std::string &output_path, GraphData<T>& graph_data, bool print=false);
};


template<typename T>
void LoadSaveGraphDatasets::LoadTUDortmundGraphData(const std::string &path, const std::string &dbName, GraphData<T> &graphs,
                                       std::vector<int> &graphLabels,
                                       std::vector<std::vector<Label>> *graphsNodeLabels,
                                       std::vector<std::vector<Label>> *graphsEdgeLabels,
                                       std::vector<std::vector<double>> *graphsNodeAttributes,
                                       std::vector<std::vector<double>> *graphsEdgeAttributes) {
    static_assert(std::is_base_of_v<GraphStruct, T>, "T must derive from GraphStruct");
    graphs.SetName(dbName);
    EDGES edges;
    std::vector<int> graphIndicator;
    std::vector<Label> nodeLabels;
    std::vector<Label> edgeLabels;
    std::vector<double> nodeAttributes;
    std::vector<double> edgeAttributes;
    bool isNodeLabels = false;
    bool isEdgeLabels = false;
    bool isNodeAttributes = false;
    bool isEdgeAttributes = false;
    //read edge data
    std::ifstream edge_file(path + dbName + "/" + dbName + "_A.txt", std::ios_base::out);
    std::string edge_string;
    std::string number_string;
    std::string str;
    if (edge_file.is_open()){
        while (std::getline(edge_file, str)){
            std::stringstream edge_stream(str);
            std::pair<INDEX, INDEX> edge;
            bool first = true;
            while (std::getline(edge_stream, number_string, ',')) {
                if (first){
                    first = false;
                    edge.first = std::stoi(number_string);
                }
                else{
                    edge.second = std::stoi(number_string);
                }
            }
            edges.push_back(edge);
        }
        edge_file.close();
    }

    //read _graph indicator data
    std::ifstream graph_indicator_file(path + dbName + "/" + dbName + "_graph_indicator.txt", std::ios_base::out);

    if (graph_indicator_file.is_open()){
        while (std::getline(graph_indicator_file, str)){
            graphIndicator.emplace_back(std::stoi(str));
        }
        graph_indicator_file.close();
    }

    //read _graph labels
    std::ifstream graph_label_file(path + dbName + "/" + dbName + "_graph_labels.txt", std::ios_base::out);
    if (graph_label_file.is_open()){
        while (std::getline(graph_label_file, str)){
            graphLabels.emplace_back(std::stoi(str));
        }
        graph_label_file.close();
    }

    //read node labels
    if(std::filesystem::is_regular_file(path + dbName + "/" + dbName + "_node_labels.txt")) {
        isNodeLabels = true;
        std::ifstream node_label_file(path + dbName + "/" + dbName + "_node_labels.txt", std::ios_base::out);
        if (node_label_file.is_open()) {
            while (std::getline(node_label_file, str)) {
                nodeLabels.emplace_back(std::stoi(str));
            }
            node_label_file.close();
        }
    }
    //read edge labels
    if(std::filesystem::is_regular_file(path + dbName + "/" + dbName + "_edge_labels.txt")) {
        isEdgeLabels = true;
        std::ifstream edge_label_file(path + dbName + "/" + dbName + "_edge_labels.txt", std::ios_base::out);
        if (edge_label_file.is_open()) {
            while (std::getline(edge_label_file, str)) {
                edgeLabels.emplace_back(std::stoi(str));
            }
        }
    }

    //read node attributes
    if(std::filesystem::is_regular_file(path + dbName + "/" + dbName + "_node_attributes.txt")) {
        isNodeAttributes = true;
        std::ifstream node_attribute_file(path + dbName + "/" + dbName + "_node_attributes.txt", std::ios_base::out);
        if (node_attribute_file.is_open()) {
            while (std::getline(node_attribute_file, str)) {
                nodeAttributes.emplace_back(std::stoi(str));
            }
            node_attribute_file.close();
        }
    }
    //read edge attributes
    if(std::filesystem::is_regular_file(path + dbName + "/" + dbName + "_edge_attributes.txt")) {
        isEdgeAttributes = true;
        std::ifstream edge_attribute_file(path + dbName + "/" + dbName + "_edge_attributes.txt", std::ios_base::out);
        if (edge_attribute_file.is_open()) {
            while (std::getline(edge_attribute_file, str)) {
                edgeAttributes.emplace_back(std::stoi(str));
            }
        }
    }

    //Get all _graph nodes
    std::unordered_map<INDEX, std::pair<INDEX, INDEX>> idMap;
    INDEX nodeCounter = 0;
    INDEX graphNodeCounter = 0;
    for (auto indicator : graphIndicator) {
        ++nodeCounter;
        ++graphNodeCounter;
        if (indicator > graphs.size()){
            graphNodeCounter = 1;
            std::string graph_name = dbName + "_" + std::to_string(indicator-1);
            graphs.add(T(graph_name, 0));
        }
        graphs[indicator-1].AddNodes(1);

        //get node labels
        if (isNodeLabels && graphsNodeLabels != nullptr) {
            if (graphsNodeLabels->size() < indicator) {
                graphsNodeLabels->resize(indicator);
            }
            (*graphsNodeLabels)[indicator - 1].emplace_back(nodeLabels[nodeCounter-1]);
        }

        //get node attributes
        if (isNodeAttributes && graphsNodeAttributes != nullptr){
            if (graphsNodeAttributes->size() < indicator) {
                graphsNodeAttributes->resize(indicator);
            }
            (*graphsNodeAttributes)[indicator - 1].emplace_back(nodeAttributes[nodeCounter-1]);
        }

        idMap.insert({nodeCounter,{graphs.size() - 1, graphNodeCounter - 1}});
    }

    // Add edge labels
    if (isEdgeLabels && graphsEdgeLabels != nullptr){
        graphsEdgeLabels->clear();
        for (auto const & graph : graphs.graphData) {
            graphsEdgeLabels->emplace_back(std::vector<Label>());
        }
    }
    if (isEdgeAttributes && graphsEdgeAttributes != nullptr){
        graphsEdgeAttributes->clear();
        for (auto const & graph : graphs.graphData) {
            graphsEdgeAttributes->emplace_back(std::vector<double>());
        }
    }

    //Get all _graph edges
    INDEX edgeCounter = 0;
    for (auto edge : edges) {
        std::vector<double> edge_data;
        INDEX graphId = idMap[edge.first].first;
        if (graphId != idMap[edge.second].first){
            throw std::logic_error("File format is wrong!");
        }
        else{
            if (isEdgeLabels && graphsEdgeLabels != nullptr) {
                //std::cout << graphId << ":"  << graphsEdgeLabels->size() << std::endl;
                //std::cout << edgeCounter << ":"  << edgeLabels.size() << std::endl;
                (*graphsEdgeLabels)[graphId].emplace_back(edgeLabels[edgeCounter]);
                edge_data.emplace_back(edgeLabels[edgeCounter]);
            }
            if (isEdgeAttributes && graphsEdgeAttributes != nullptr) {
                (*graphsEdgeAttributes)[graphId].emplace_back(edgeAttributes[edgeCounter]);
                edge_data.emplace_back(edgeAttributes[edgeCounter]);
            }
            graphs[graphId].AddEdge(idMap[edge.first].second, idMap[edge.second].second, edge_data);
        }
        ++edgeCounter;
    }

    // check whether T can be casted to UDataGraph
    if (std::is_same_v<T, UDataGraph>) {
        int graph_counter = 0;
        for (auto &graph : graphs.graphData) {
            // cast
            auto & data_graph = dynamic_cast<UDataGraph&>(graph);
            std::vector<std::string> nodeFeatureNames = {"label"};
            if (isNodeAttributes) {
                int attr_counter = 1;
                for (auto const & node : graphsNodeAttributes[0]) {
                    std::string attr_name = "attr_" + std::to_string(attr_counter);
                    nodeFeatureNames.emplace_back(attr_name);
                    ++attr_counter;
                }
            }
            std::vector<std::string> edgeFeatureNames;
            if (isEdgeLabels) {
                edgeFeatureNames.emplace_back("label");
            }
            if (isEdgeAttributes) {
                int attr_counter = 1;
                for (auto const & edge : graphsEdgeAttributes[0]) {
                    std::string attr_name = "attr_" + std::to_string(attr_counter);
                    edgeFeatureNames.emplace_back(attr_name);
                    ++attr_counter;
                }
            }
            data_graph.Init(data_graph.GetName(),
                data_graph.nodes(), data_graph.edges(),
                nodeFeatureNames.size(),
                edgeFeatureNames.size(),
                nodeFeatureNames,
                edgeFeatureNames);

            std::vector<double> node_data;
            for (INDEX i = 0; i < data_graph.nodes(); ++i) {
                node_data.clear();
                // Initialize graph node labels
                if (isNodeLabels && graphsNodeLabels != nullptr) {
                    data_graph.ReadNodeFeatures((*graphsNodeLabels)[graph_counter][i], i, "label");
                } else {
                    data_graph.ReadNodeFeatures(0, i, "label");
                }
                // add node attributes
                if (isNodeAttributes && graphsNodeAttributes != nullptr) {
                    int attr_counter = 1;
                    for (auto const & attr : (*graphsNodeAttributes)[idMap[i + 1].first]) {
                        data_graph.ReadNodeFeatures(attr, i, "attr_" + std::to_string(attr_counter));
                        ++attr_counter;
                    }
                }
            }
            std::cout << "Graph " << data_graph.GetName() << " has " << data_graph.nodes() << " nodes and " << data_graph.edges() << " edges." << std::endl;
            std::cout << "Number of Labels: " << (*graphsNodeLabels)[graph_counter].size() << std::endl;
            data_graph.SetLabels(&(*graphsNodeLabels)[graph_counter]);
            std::cout << "Set node labels for graph " << data_graph.GetName() << std::endl;
            ++graph_counter;
        }
    }
}

inline bool LoadSaveGraphDatasets::PreprocessTUDortmundGraphData(const std::string &dataset_name, const std::string &input_path, const std::string &output_path) {
     // check whether folder with name dataset_name is not yet existing
    if (!std::filesystem::exists(input_path + dataset_name + "/") && !std::filesystem::is_directory(output_path + dataset_name))
    {
        std::cout << "Folder " << input_path << " does not exist" << std::endl;
        // Please create it and download the data dataset_name from link
        std::cout << "Please download the dataset from https://chrsmrrs.github.io/datasets/docs/datasets/" << std::endl;
        return false;
    }
    // if output_path except for ProcessedGraphs exists then create ProcessedGraphs folder
    if (!std::filesystem::exists(output_path)) {
        std::filesystem::create_directory(output_path);
    }



    // Check whether Processed graphs already exists
    if (std::filesystem::exists(output_path + dataset_name + ".bgf")) {
        std::cout << "Graph " << dataset_name << " already exists" << std::endl;
        return true;
    }

    UDataGraph graph;
    GraphData<UDataGraph> graphs;
    //graphs.add(example_graph());
    std::vector<int> graphLabels;
    std::vector<std::vector<INDEX>> graphNodeLabels;
    std::vector<std::vector<double>> graphNodeAttributes;
    std::vector<std::vector<double>> graphEdgeAttributes;
    std::vector<std::vector<INDEX>> graphEdgeLabels;
    LoadSaveGraphDatasets::LoadTUDortmundGraphData(input_path, dataset_name, graphs, graphLabels, &graphNodeLabels, &graphEdgeLabels, &graphNodeAttributes, &graphEdgeAttributes);

    SaveParams params = {
        output_path,
        dataset_name,
        GraphFormat::BGF,
        true,
    };
    // Save the graph as bgfs format
    graphs.Save(params);
    return true;
}

template <typename T>
void LoadSaveGraphDatasets::LoadPreprocessedTUDortmundGraphData(const std::string& dataset_name , const std::string &output_path, GraphData<T>& graph_data, bool print) {
    // base class should be GraphStruct
    static_assert(std::is_base_of_v<GraphStruct, T>, "T must derive from GraphStruct");
    // Load the graph from the bgfs format
    std::string graph_path = output_path + dataset_name + ".bgf";

    if (!std::filesystem::exists(graph_path)) {
        std::cout << "Graph " << dataset_name << " does not exist" << std::endl;
    }
    else {
        graph_data.Load(graph_path);
        // print the loaded graphs
        if (print) {
            for ( auto &x : graph_data.graphData) {
                std::cout << x << std::endl;
            }
        }
        std::cout << "Successfully loaded the " << dataset_name << " graphs from TUDataset" << std::endl;
    }
    graph_data.SetName(dataset_name);
}


#endif //LOADSAVEGRAPHDATASETS_H
