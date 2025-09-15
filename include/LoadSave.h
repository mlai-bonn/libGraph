//
// Created by Florian on 16.03.2021.
//

#ifndef LOADSAVE_H
#define LOADSAVE_H

#include <fstream>
#include <string>
#include <set>
#include <unordered_map>
#include "typedefs.h"
#include "GraphDataStructures/GraphStructs.h"


class LoadSave {
public:
    static bool LoadLabels(const std::string& Path, Labels &labels);
    static void LoadLabelsFromPath(const std::string& Path, std::vector<Labels> &Labels, std::unordered_map<size_t, std::string> &LabelNames, std::set<int>* graphSizes = nullptr, int patternNum = -1);


    /// Converts all graphs in a directory using the same extension to the given format
    /// \param path
    /// \param Format
    template<typename T>
    static void Convert(const std::string & path, GraphFormat Format = GraphFormat::BGF, const std::string & extension = ".txt");

    static bool ReadBGF(const std::string &extension, std::ifstream &In, int &saveVersion, int &graphNumber,
                        std::vector<std::string> &graphsNames, std::vector<GraphType> &graphsTypes,
                        std::vector<INDEX> &graphsSizes, std::vector<std::vector<std::string>> &graphsNodeFeatureNames,
                        std::vector<INDEX> &graphsEdges, std::vector<std::vector<std::string>> &graphsEdgeFeatureNames);

    template<typename  T>
    static void LoadGraphsFromPath(const std::string& graphPath, const std::string& labelPath, std::vector<T> &graphs, const std::string & searchName = "", const std::string & extension = "", bool sort = true, std::set<int>* graphSizes = nullptr, int patternNum = -1);
    template<typename T>
    static void Convert(const std::string & graphPath, const std::string& Name = "", GraphFormat Format = GraphFormat::BGF, bool relabeling = true, bool withLabels = false, const std::string& labelPath = "", bool Labeled = false);

};

inline bool LoadSave::LoadLabels(const std::string& Path, Labels &labels) {
    labels.clear();
    if (std::filesystem::path(Path).extension() == ".labels") {
        try {
            std::ifstream infile(Path);
            NodeId id;
            Label label;
            while (infile >> id >> label){
                labels.push_back(label);
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
inline void LoadSave::Convert(const std::string &graphPath, const std::string &Name,
                                 GraphFormat Format, bool relabeling, bool withLabels, const std::string &labelPath, bool Labeled) {
    T graph = T(graphPath, relabeling, withLabels, labelPath);
    graph.Save({"", Name, Format, Labeled});
}

template<typename T>
inline void LoadSave::Convert(const std::string &path, GraphFormat Format, const std::string &extension) {
    std::vector<std::string> files;
    for (const auto &entry: std::filesystem::directory_iterator(path)) {
        if (entry.is_regular_file()) {
            files.emplace_back(entry.path().string());
        }
    }
    for(const auto& f_path : files){
        if (std::filesystem::path(f_path).extension() == extension){
            T graph = T(f_path, true, false, "");
            graph.Save({"", "", Format, false});
        }
    }

}
template<typename T>
inline void LoadSave::LoadGraphsFromPath(const std::string& graphPath, const std::string& labelPath, std::vector<T> &graphs, const std::string & searchName, const std::string & extension, const bool sort, std::set<int>* graphSizes, int patternNum)
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
        std::sort(graphs.begin(), graphs.end(), T::sort_by_name());
    }

}

inline bool LoadSave::ReadBGF(const std::string& extension,std::ifstream& In, int &saveVersion,  int &graphNumber, std::vector<std::string> &graphsNames,
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


#endif //LOADSAVE_H
