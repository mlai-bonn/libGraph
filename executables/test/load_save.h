//
// Created by florian on 23.01.23.
//

#ifndef TESTGRAPHLIB_UNITTEST_LOAD_SAVE_H
#define TESTGRAPHLIB_UNITTEST_LOAD_SAVE_H


#include <string>

#include "load_save.h"
#include "../../include/GraphDataStructures/GraphBase.h"
#include "../../include/io/StaticFunctions.h"

bool TestLoadGraphsFromPath(const std::string& graph_path, const std::string& label_path, const std::string& extension) {

    std::vector<GraphStruct> graphs;
    GraphStruct::LoadGraphsFromPath(graph_path, label_path, graphs, extension);
    return graphs.size() > 0;
}

void ConvertDimacs(){
    for (const auto& e1: std::filesystem::directory_iterator("../../../../GraphData/Hops/Automorphisms/")) {
        if (is_directory(std::filesystem::path(e1))) {
            GraphData graphData = GraphData<GraphStruct>();
            std::string stem = std::filesystem::path(e1).stem().string();
            for (const auto &e2: std::filesystem::recursive_directory_iterator(e1.path())) {
                std::string path = e2.path().string();
                std::cout << path << std::endl;
                GraphStruct g = GraphStruct(path, true, false, "", "dimacs");
                graphData.add(g);
            }
            graphData.Save({"../../../../GraphData/Hops/Automorphisms/", stem, GraphFormat::BGFS});
            //GraphData g = GraphData<GraphStruct>("../../../../GraphData/Hops/Automorphisms/cfi-rigid-d3.bgfs");
            //int num = g.graphData.size();
        }
    }
}

void LoadAids(){
    GraphData g = GraphData<GraphStruct>(GraphFormat::AIDS, "../../../../GraphData/Hops/patterns/labeled/tree_v5e1_levels/tree_v5e1_size2.graph");
    return;
}

void LoadSpeedTest(){
    auto start = std::chrono::high_resolution_clock::now();
    GraphStruct graphStruct = GraphStruct("../../../GraphData/Hops/com-lj.ungraph.bgfs");
    auto end = std::chrono::high_resolution_clock::now();
    INDEX size = graphStruct.nodes();
    INDEX edges = graphStruct.edges();
    std::cout << "Loaded " << size << " nodes " << " and " << edges << " edges " << " in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000.0 << "s" << std::endl;
}

bool LoadSaveTest(){
    DDataGraph dataGraph = DDataGraph("../test/TestData/DGraph.edges", true, true,"../test/TestData/DGraph.data");
    dataGraph.Save({"../test/TestData/", "DGraph"});
    DDataGraph loadGraph = DDataGraph("../test/TestData/DGraph.bgf");
    auto graphData = GraphData<DDataGraph>("../test/TestData/DGraph.bgf");
    bool x = (dataGraph == loadGraph);
    x = x && (dataGraph == graphData.graphData[0]);

    graphData.Save({"../test/TestData/", "XGraph"});
    loadGraph = DDataGraph("../test/TestData/DGraph.bgf");
    x = x && (dataGraph == loadGraph);

    DGraphStruct dgraphStruct = DGraphStruct("../test/TestData/DGraph.bgf");
    auto dgraphData = GraphData<DGraphStruct>("../test/TestData/DGraph.bgf");
    x = x && (dgraphStruct == dgraphData.graphData[0]);
    dgraphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGF});
    dgraphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGFS});

    DGraphStruct dbgf = DGraphStruct("../test/TestData/DGraph.bgf");
    DGraphStruct dbgfs = DGraphStruct("../test/TestData/DGraph.bgfs");

    x = x && (dgraphStruct == dbgf);
    x = x && (dbgf == dbgfs);

    dgraphData.Save({"../test/TestData/", "DGraph", GraphFormat::BGF});

    dbgf = DGraphStruct("../test/TestData/DGraph.bgf");

    x = x && (dgraphData.graphData[0] == dbgf);

    GraphStruct graphStruct = GraphStruct("../test/TestData/DGraph.bgf");
    graphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGF});
    graphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGFS});

    GraphStruct bgf = GraphStruct("../test/TestData/DGraph.bgf");
    GraphStruct bgfs = GraphStruct("../test/TestData/DGraph.bgfs");


    x = x && (graphStruct == bgf);
    x = x && (bgf == bgfs);

    dgraphStruct = DGraphStruct("../test/TestData/DGraph.bgf");
    dgraphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGF});
    dgraphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGFS});

    dbgf = DGraphStruct("../test/TestData/DGraph.bgf");
    dbgfs = DGraphStruct("../test/TestData/DGraph.bgfs");

    x = x && (dgraphStruct == dbgf);
    x = x && (dbgf == dbgfs);

    dataGraph = DDataGraph("../test/TestData/DGraph.edges", true, true,"../test/TestData/DGraph.data");
    dataGraph.Save({"../test/TestData/", "DGraph"});
    loadGraph = DDataGraph("../test/TestData/DGraph.bgf");
    x = x || (dataGraph == loadGraph);

    graphStruct = GraphStruct("../test/TestData/DGraph.bgf");
    graphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGF});
    graphStruct.Save({"../test/TestData/", "DGraph", GraphFormat::BGFS});

    bgf = GraphStruct("../test/TestData/DGraph.bgf");
    bgfs = GraphStruct("../test/TestData/DGraph.bgfs");


    x = x && (graphStruct == bgf);
    x = x && (bgf == bgfs);

    return x;
}

bool LoadSaveTest2(){
    for (int i = 3; i < 5; ++i) {
        GraphData graphData = GraphData<GraphStruct>();
        graphData.Load("../../../GraphData/Hops/patterns/size" + std::to_string(i) + ".bgf");
        graphData.add("../../../GraphData/Hops/patterns/size" + std::to_string(i) + "/", "", "", ".txt");
        graphData.Save({"../../../GraphData/Hops/patterns/", "size" + std::to_string(i)});
        graphData.Load("../../../GraphData/Hops/patterns/size" + std::to_string(i) + "/size" + std::to_string(i) + ".bgf");
        int x = graphData.size();
    }
}


bool LoadCSV(){
    std::vector<std::vector<std::string>> out;
    StaticFunctionsLib::load_csv("../../../GraphData/Hops/patterns/automorphisms_3.csv", out);
    std::vector<int> automorphisms;
    int counter = 0;
    for (auto &x : out) {
        if (counter > 0) {
            automorphisms.emplace_back(std::stoi(x[3]));
        }
        ++counter;
    }
    return true;
}

std::map<std::string, GraphStruct> LoadTXTGraphs(){
    std::vector<std::string> paths = StaticFunctionsLib::directory_paths("../../../../GraphData/data_ravkic/graphs/", {"DBLP", "FACEBOOK", "WEBKB", "YEAST"});
    GraphStruct graphStruct;
    std::map<std::string, GraphStruct> graphStructMap;
    for (auto &path : paths) {
            // iterate over path and get all txt files in it
            std::vector<std::string> files = StaticFunctionsLib::file_paths(path, {}, {".txt"});
            // Add path to graphMap keys
            for (auto &file : files){
                // load graph from edge list
                graphStructMap[path] = GraphStruct(file);
            }
    }
    return graphStructMap;
}


std::map<std::string, std::map<int, GraphData<GraphStruct>>> LoadTXTPatterns(){
    std::vector<std::string> paths = StaticFunctionsLib::directory_paths("../../../../GraphData/data_ravkic/patterns/", {"DBLP", "FACEBOOK", "WEBKB", "YEAST"});
    GraphData<GraphStruct> graphData;
    std::map<std::string, std::map<int, GraphData<GraphStruct>>> graphDataMap;
    for (auto &path : paths) {
            // iterate over path and get all txt files in it
            graphDataMap.emplace(path, std::map<int, GraphData<GraphStruct>>());
            std::vector<std::string> files = StaticFunctionsLib::file_paths(path, {}, {".txt"});
            // Add path to graphMap keys
            for (auto &file : files){
                // load graph from edgelist
                GraphStruct graphStruct = GraphStruct(file);
                int node_number = graphStruct.nodes();
                // check if node number is already in graphMap
                if (graphDataMap[path].find(node_number) == graphDataMap[path].end()){
                    // if not add it
                    graphDataMap[path].emplace(node_number, GraphData<GraphStruct>());

                }
                graphDataMap[path][node_number].add(graphStruct);
                // the id is the size of the curren graphMap entry
                int id = graphDataMap[path][node_number].size();
                // set name path, node number, edge number and counter separated by _
                graphDataMap[path][node_number].graphData.back().SetName(std::to_string(id));
            }
        }
    for (auto &x : graphDataMap) {
        // print path
        for (auto &y : x.second) {
            std::cout << "Graph:" << std::endl;
            std::cout << x.first << std::endl;
            std::cout << "Node number:" << std::endl;
            std::cout << y.first << std::endl;
            for (auto & graph : y.second.graphData){
                // print path : number of nodes : number of edges
                std::cout << graph.GetName() << " : ";
                std::cout << graph.nodes() << " : ";
                std::cout << graph.edges() << " : ";
                std::cout << std::endl;
            }
        }
    }
    return graphDataMap;
}

#endif //TESTGRAPHLIB_UNITTEST_LOAD_SAVE_H
