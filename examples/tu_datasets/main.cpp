#include <iostream>
#include <cstring>
#include <chrono>
#include <filesystem>
#include <vector>
#include <LoadSave.h>
#include <GraphFunctions.h>
#include <omp.h>

// hello world

int main() {
    std::string dataset_name = "MUTAG";
    // Load the MUTAG graph
    std::string input_path = "../Data/";
    std::string output_path = "../Data/ProcessedGraphs/";

    // check wheter folder "../Data/" exsits
    if (!std::filesystem::exists(input_path)) {
        std::cout << "Folder " << input_path << " does not exist" << std::endl;
        // Please create it and download the data dataset_name from link
        std::cout << "Please download the dataset from https://chrsmrrs.github.io/datasets/docs/datasets/" << std::endl;
        return 1;
    }
    else {
        // if output_path except for ProcessedGraphs exists then create ProcessedGraphs folder
        if (!std::filesystem::exists(output_path)) {
            std::filesystem::create_directory(output_path);
        }
    }




    GraphStruct graphStruct;
    GraphData<GraphStruct> graphs;
    //graphs.add(example_graph());
    std::vector<int> graphLabels;
    std::vector<std::vector<int>> graphNodeLabels;
    std::vector<std::vector<int>> graphNodeAttributes;
    std::vector<std::vector<int>> graphEdgeAttributes;
    std::vector<std::vector<int>> graphEdgeLabels;
    LoadSave::LoadTUDortmundGraphData(input_path, dataset_name, graphs, graphLabels, &graphNodeLabels, &graphEdgeLabels, &graphNodeAttributes, &graphEdgeAttributes);

    SaveParams params = {
        output_path,
        dataset_name,
        GraphFormat::BGFS,
        true,
    };
    // Add information to graphs
    int counter = 0;
    for ( auto &x : graphs.graphData) {
        auto labels = Labels(graphNodeLabels[counter].begin(), graphNodeLabels[counter].end());
        std::string graph_name = dataset_name + "_" + std::to_string(counter);
        x.SetName(graph_name);
        x.set_labels(&labels);
        ++counter;
        // print counter
        std::cout << "Graph:" << counter << "/" << graphs.graphData.size() << std::endl;
        std::cout << x.nodes() << std::endl;
    }
    // Save the graph as bgfs format
    graphs.Save(params);

    // Load the graph from the bgfs format

    GraphData<GraphStruct> loadedGraphs;
    std::string graph_path = output_path + dataset_name + ".bgfs";
    loadedGraphs.Load(graph_path);

    // print the loaded graphs
    for ( auto &x : loadedGraphs.graphData) {
        std::cout << x << std::endl;
    }

    std::cout << "Successfully loaded the graphs from TUDataset" << dataset_name << std::endl;
    return 0;
}