//
// Created by Florian on 16.03.2021.
//

#ifndef HOPS_LOADSAVE_H
#define HOPS_LOADSAVE_H

#include <string>
#include <set>
#include <unordered_map>
#include "typedefs.h"
#include "DataClasses.h"


class LoadSave {
public:
    static bool LoadLabels(const std::string& Path, Labels &Labels);

    static void LoadLabelsFromPath(const std::string& Path, std::vector<Labels> &Labels, std::unordered_map<size_t, std::string> &LabelNames, std::set<int>* graphSizes = nullptr, int patternNum = -1);
    template <typename T>
    static void LoadTUDortmundGraphData(const std::string& path, const std::string& dbName, GraphData<T>& graphs, std::vector<int>& graphLabels, std::vector<std::vector<int>>* graphsEdgeLabels= nullptr, std::vector<std::vector<int>>* graphsNodeAttributes= nullptr, std::vector<std::vector<int>>* graphsEdgeAttributes= nullptr);

};

inline bool LoadSave::LoadLabels(const std::string& Path, Labels &LabelVector) {
    LabelVector.clear();
    if (std::filesystem::path(Path).extension() == ".labels") {
        try {
            std::ifstream infile(Path);
            NodeId id;
            Label label;
            while (infile >> id >> label){
                LabelVector.push_back(label);
            }
        }
        catch (...) {
            return false;
        }
        return true;
    }
    return false;
}

void
inline LoadSave::LoadLabelsFromPath(const std::string& Path, std::vector<Labels> &Labels, std::unordered_map<size_t, std::string> &LabelNames, std::set<int>* graphSizes, int patternNum) {
    std::vector<Label> LabelVector;
    std::unordered_map<INDEX, INDEX> sizesNumMap;
    if (graphSizes != nullptr) {
        for (INDEX size : *graphSizes) {
            sizesNumMap.insert({size, 0});
        }
    }

    for (const auto &entry : std::filesystem::directory_iterator(Path)) {
        if (LoadLabels(entry.path().string(), LabelVector)) {
            LabelNames.insert(std::pair<size_t, std::string>(Labels.size(), entry.path().stem().string()));
            if (graphSizes == nullptr || graphSizes->find(LabelVector.size()) != graphSizes->end()) {
                if (patternNum == -1 || graphSizes == nullptr || sizesNumMap[LabelVector.size()] < patternNum) {
                    Labels.push_back(LabelVector);
                    if (graphSizes != nullptr) {
                        ++sizesNumMap[LabelVector.size()];
                    }
                }
            }

        }
    }
}

template<typename T>
void LoadSave::LoadTUDortmundGraphData(const std::string &path, const std::string &dbName, GraphData<T> &graphs,
                                       std::vector<int> &graphLabels, std::vector<std::vector<int>> *graphsEdgeLabels,
                                       std::vector<std::vector<int>> *graphsNodeAttributes,
                                       std::vector<std::vector<int>> *graphsEdgeAttributes) {

    EDGES edges;
    std::vector<int> graphIndicator;
    std::vector<int> nodeLabels;
    std::vector<int> edgeLabels;
    std::vector<int> nodeAttributes;
    std::vector<int> edgeAttributes;
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

    //read graph indicator data
    std::ifstream graph_indicator_file(path + dbName + "/" + dbName + "_graph_indicator.txt", std::ios_base::out);

    if (graph_indicator_file.is_open()){
        while (std::getline(graph_indicator_file, str)){
            graphIndicator.emplace_back(std::stoi(str));
        }
        graph_indicator_file.close();
    }

    //read graph labels
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
    //read node attributes
    if(std::filesystem::is_regular_file(path + dbName + "/" + dbName + "_edge_attributes.txt")) {
        isEdgeAttributes = true;
        std::ifstream edge_attribute_file(path + dbName + "/" + dbName + "_edge_attributes.txt", std::ios_base::out);
        if (edge_attribute_file.is_open()) {
            while (std::getline(edge_attribute_file, str)) {
                edgeAttributes.emplace_back(std::stoi(str));
            }
        }
    }

    //Get all graph nodes
    std::unordered_map<INDEX, std::pair<INDEX, INDEX>> idMap;
    INDEX nodeCounter = 0;
    INDEX graphNodeCounter = 0;
    for (auto indicator : graphIndicator) {
        ++nodeCounter;
        ++graphNodeCounter;
        Labels labels;
        if (indicator > graphs.size()){
            graphNodeCounter = 1;
            graphs.add(GraphStruct());
            if (!labels.empty()){
                if (isNodeLabels){
                    graphs.graphData.back().set_labels(&labels);
                }
            }
        }
        graphs[indicator-1].add_node();

        //get node labels
        if (isNodeAttributes && graphsNodeAttributes != nullptr){
            if (graphsNodeAttributes->size() < indicator) {
                graphsNodeAttributes->resize(indicator);
            }
            (*graphsNodeAttributes)[indicator - 1].emplace_back(nodeAttributes[nodeCounter-1]);
        }
        idMap.insert({nodeCounter,{graphs.size() - 1, graphNodeCounter - 1}});
    }


    if (isEdgeLabels && graphsEdgeLabels != nullptr){
        graphsEdgeLabels->clear();
        for (auto const & graph : graphs.graphData) {
            graphsEdgeLabels->emplace_back(std::vector<int>());
        }
    }
    if (isEdgeAttributes && graphsEdgeAttributes != nullptr){
        graphsEdgeAttributes->clear();
        for (auto const & graph : graphs.graphData) {
            graphsEdgeAttributes->emplace_back(std::vector<int>());
        }
    }

    //Get all graph edges
    INDEX edgeCounter = 0;
    for (auto edge : edges) {
        INDEX graphId = idMap[edge.first].first;
        if (graphId != idMap[edge.second].first){
            throw std::logic_error("File format is wrong!");
        }
        else{
            graphs[graphId].add_edge(idMap[edge.first].second, idMap[edge.second].second);
            if (isEdgeLabels && graphsEdgeLabels != nullptr) {
                //std::cout << graphId << ":"  << graphsEdgeLabels->size() << std::endl;
                //std::cout << edgeCounter << ":"  << edgeLabels.size() << std::endl;
                (*graphsEdgeLabels)[graphId].emplace_back(edgeLabels[edgeCounter]);
            }
            if (isEdgeAttributes && graphsEdgeAttributes != nullptr) {
                (*graphsEdgeAttributes)[graphId].emplace_back(edgeAttributes[edgeCounter]);
            }
        }
        ++edgeCounter;
    }
}




#endif //HOPS_LOADSAVE_H
