//
// Created by florian on 06.11.23.
//


#include <string>
#include <filesystem>
#include <cstring>
#include <omp.h>
#include "libGraph.h"


// load all graphs from a path and convert them into bgfs format
void convert_graphs(const std::string& path, const int threads, const std::string& extension, const std::string& format, std::string& output_path){
    // recursively iterate over all files in the given path
    for (const auto & entry : std::filesystem::recursive_directory_iterator(path)){
        auto ext = entry.path().extension().string();
        // if extension is found in the file name
        if (ext.find(extension) == std::string::npos){
            continue;
        }
        // check if the file is a graphml file
        if (entry.path().extension() == ".graphml"){
            // load the graph
            GraphData<GraphStruct> graphs(entry.path().string());
            graphs.Save({.graphPath = output_path, .Name = entry.path().filename().string(), .Format = GraphFormat::BGFS});
        }
        GraphData<GraphStruct> graphs(entry.path().string());

        // set graph format using format
        // check if format contains bgfs
        if (format.find("bgfs") != std::string::npos){
            // extension is bgfs
            graphs.Save({.graphPath = output_path, .Name = entry.path().filename().string(), .Format = GraphFormat::BGFS});
        }
        // check if format contains edges
        if (format.find("edges") != std::string::npos){
            if (output_path.empty())
            {
                output_path = entry.path().parent_path().string() + "/";
            }
            graphs.Save({.graphPath = output_path, .Name = entry.path().filename().string(), .Format = GraphFormat::EDGES});
        }
    }
}

// main
int main(int argc, char* argv[]){
    // -p: path to the graphs
    // -j: number of threads
    // -e: file extension of the graphs
    // -f: output format
    // -o: output path

    std::string path;
    int threads = omp_get_max_threads();
    std::string extension = ".edges";
    std::string format = "bgfs";
    std::string output_path;

    for (int i = 0; i < argc; ++i) {
        // help
        if (strcmp(argv[i], "-h") == 0){
            std::cout << "Converts all graphs in the given path to the given output format." << std::endl;
            std::cout << "Arguments:" << std::endl;
            std::cout << "-p: path to the graphs" << std::endl;
            std::cout << "-j: number of threads (default: maximum number of threads)" << std::endl;
            std::cout << "-e: file extension of the graphs (default: .edges)" << std::endl;
            std::cout << "-f: output format (default: bgfs)" << std::endl;
            std::cout << "-o: output path (default: empty, i.e., the same as the input path)" << std::endl;
            return 0;
        }
        // read the arguments
        if (strcmp(argv[i], "-p") == 0 && i + 1 < argc){
            path = argv[i + 1];
        }
        if (strcmp(argv[i], "-j") == 0 && i + 1 < argc){
            threads = std::stoi(argv[i + 1]);
        }
        if (strcmp(argv[i], "-e") == 0 && i + 1 < argc){
            extension = argv[i + 1];
        }
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc){
            format = argv[i + 1];
        }
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc){
            output_path = argv[i + 1];
        }
    }
    convert_graphs(path, threads, extension, format, output_path);
    return 0;
}