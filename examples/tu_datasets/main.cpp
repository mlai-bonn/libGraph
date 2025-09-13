#include <iostream>
#include <chrono>
#include <filesystem>
#include <vector>
#include <LoadSave.h>

// hello world

int main() {
    std::string dataset_name = "MUTAG";
    // Load the MUTAG graph
    const std::string data_path = "../Data/";
    const std::string input_path = "../Data/Graphs/";
    std::string output_path = "../Data/ProcessedGraphs/";

    // check whether folder "../Data/" exsits
    if (!std::filesystem::exists(input_path + dataset_name))
    {
        std::cout << "Input path " << input_path + dataset_name << " does not exist." << std::endl;
        // Please create it and download the data dataset_name from link
        std::cout << "Please download the dataset from https://chrsmrrs.github.io/datasets/docs/datasets/" << std::endl;
        return 1;
    }
    // if output_path except for ProcessedGraphs exists then create ProcessedGraphs folder
    if (std::filesystem::exists(data_path) && !std::filesystem::exists(output_path)) {
        std::filesystem::create_directory(output_path);
    }

    LoadSave::PreprocessTUDortmundGraphData(dataset_name, input_path, output_path);
    GraphData<UDataGraph> loadedGraphs;
    LoadSave::LoadPreprocessedTUDortmundGraphData(dataset_name, output_path, loadedGraphs);

    // print the loaded graphs
    for ( auto &x : loadedGraphs.graphData) {
        std::cout << x << std::endl;
    }

    std::cout << "Successfully loaded the " << dataset_name << " graphs from TUDataset" << std::endl;
    return 0;
}