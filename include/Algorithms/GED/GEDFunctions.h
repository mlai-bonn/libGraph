//
// Created by florian on 29.08.25.
//

#ifndef GEDEXAMPLE_GEDFUNCTIONS_H
#define GEDEXAMPLE_GEDFUNCTIONS_H
#include <map>

#include "GEDEvaluation.h"
#include "typedefs.h"

// for approximated GED computation we use the gedlib library see https://github.com/dbblumenthal/gedlib
// the following functions are only to embedd the library results in our code base, i.e., the evaluation of the result

// main function to run the examples

inline std::vector<GraphStruct> CreateEditPath(const GEDEvaluation& result)
{
    std::cout << "Time: " << result.time << " seconds" << std::endl;
    // print node mapping
    //std::cout << "Node Mapping First: " << std::endl;
    //for (const auto& x : result.node_mapping.first) {
    //    std::cout << x << " ";
    //}
    //std::cout << std::endl;
    //std::cout << "Node Mapping Second: " << std::endl;
    //for (const auto& x : result.node_mapping.second) {
    //    std::cout << x << " ";
    //}
    //std::cout << std::endl;

    EditPath edit_path;
    result.get_edit_path(edit_path, 0);
    // print edit path length
    std::cout << "Edit Path Length: " << edit_path.edit_path_graphs.size() - 1 << std::endl;
    std::vector<GraphStruct> edit_path_graphs;
    for (const auto& g : edit_path.edit_path_graphs) {
        edit_path_graphs.emplace_back(g);
    }
    edit_path_graphs.back().SetName(edit_path.target_graph.GetName());
    // print all the graphs
    //for (const auto& g : edit_path_graphs.graphData) {
    //    std::cout << g << std::endl;
    //}
    return edit_path_graphs;
}

inline void CreateAllEditPaths(std::vector<GEDEvaluation> &results, GraphData<GraphStruct> &graph_data, const std::string &edit_path_output = "../Data/EditPaths/") {
    // check whether file already exists
    if (std::filesystem::exists(edit_path_output + graph_data.GetName() + "_edit_paths.bgf")) {
        std::cout << "Edit paths for " << graph_data.GetName() << " already exist." << std::endl;
        return;
    }
    GraphData<GraphStruct> all_path_graphs;
    // counter for number of computed paths
   int counter = 0;
    // time variable
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (const auto& result : results) {
        std::cout << "Computing Path between graph " << result.graph_ids.first << " and graph " << result.graph_ids.second << std::endl;
        // estimated time in minutes
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        const double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
        const double estimated_total_time = (elapsed_seconds / (counter + 1)) * results.size();
        const double estimated_time_left = estimated_total_time - elapsed_seconds;
        std::cout << "Estimated time left: " << estimated_time_left / 60 << " minutes" << std::endl;
        ++counter;
        auto edit_path_graphs = CreateEditPath(result);
        int path_counter = 0;
        for (const auto& g : edit_path_graphs) {
            g.SetName(graph_data.GetName() + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(path_counter));
            all_path_graphs.add(g);
            ++path_counter;
        }
    }

    for (int i = 0; i < graph_data.size(); ++i) {
        std::cout << graph_data[i] << std::endl;
        GraphData<GraphStruct> result;
        for (int j = i+1; j < graph_data.size(); ++j) {
            std::cout << "Computing Path between graph " << i << " and graph " << j << std::endl;
            // print percentage
            std::cout << "Progress: " << (counter * 100) / (graph_data.size() * (graph_data.size() - 1) / 2) << "%" << std::endl;
            // estimated time in minutes
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            const double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            const double estimated_total_time = (elapsed_seconds / (counter + 1)) * (graph_data.size() * (graph_data.size() - 1) / 2);
            const double estimated_time_left = estimated_total_time - elapsed_seconds;
            std::cout << "Estimated time left: " << estimated_time_left / 60 << " minutes" << std::endl;
            int path_counter = 0;
            for (auto const& result : results) {
                GraphStruct new_graph = CreateEditPath(result);
                g.SetName(graph_data.GetName() + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(path_counter));
                result.add(g);
                ++path_counter;
            }
            ++counter;
        }
        // save intermediate result in the tmp folder under datasetname_i.bgfs
        SaveParams params = {
            edit_path_output + "tmp/",
            graph_data.GetName() + "_" + std::to_string(i),
            GraphFormat::BGF,
            true,
        };
        result.Save(params);
        std::cout << "Saved intermediate result for graph " << i << std::endl;
    }
    // Merge all graphs in tmp folder
    GraphData<GraphStruct> final_result;
    for (int i = 0; i < graph_data.size(); ++i) {
        GraphData<GraphStruct> result;
        std::string file_path = edit_path_output + "tmp/" + graph_data.GetName() + "_" + std::to_string(i) + ".bgf";
        result.Load(file_path);
        for (const auto& g : result.graphData) {
            final_result.add(g);
        }
    }
    // Save final result
    SaveParams params = {
        edit_path_output,
        graph_data.GetName() + "_edit_paths",
        GraphFormat::BGF,
        true,
    };
    final_result.Save(params);
    // remove all databasename related files in tmp folder
    for (const auto& entry : std::filesystem::directory_iterator(edit_path_output + "tmp/")) {
        if (entry.path().string().find(graph_data.GetName() + "_") != std::string::npos) {
            std::filesystem::remove(entry.path());
        }
    }
}

inline void GEDResultToBinary(const std::string &target_path, GEDEvaluation &result) {
    // create output target_path if it does not exist
    if (!std::filesystem::exists(target_path)) {
        std::filesystem::create_directory(target_path);
    }
    int source_id = result.graph_ids.first;
    int target_id = result.graph_ids.second;
    const std::vector<NodeId>& i_to_j = result.node_mapping.first;
    const std::vector<NodeId>& j_to_i = result.node_mapping.second;
    // open binary file
    std::string file_path = target_path + result.graph_data_name + "_" + std::to_string(source_id) + "_" + std::to_string(target_id) + "_ged_mapping.bin";
    std::ofstream ofs(file_path, std::ios::binary | std::ios::app);
    if (!ofs) {
        std::cerr << "Error opening file for writing: " << file_path << std::endl;
        return;
    }
    // write source_id and target_id
    ofs.write(reinterpret_cast<const char *>(&source_id), sizeof(source_id));
    ofs.write(reinterpret_cast<const char *>(&target_id), sizeof(target_id));
    // write size of i_to_j
    size_t size_i_to_j = i_to_j.size();
    ofs.write(reinterpret_cast<const char *>(&size_i_to_j), sizeof(size_i_to_j));
    // write i_to_j
    ofs.write(reinterpret_cast<const char *>(i_to_j.data()), size_i_to_j * sizeof(NodeId));
    // write size of j_to_i
    size_t size_j_to_i = j_to_i.size();
    ofs.write(reinterpret_cast<const char *>(&size_j_to_i), sizeof(size_j_to_i));
    // write j_to_i
    ofs.write(reinterpret_cast<const char *>(j_to_i.data()), size_j_to_i * sizeof(NodeId));
    // write graph_database_name
    size_t name_size = result.graph_data_name.size();
    ofs.write(reinterpret_cast<const char *>(&name_size), sizeof(name_size));
    ofs.write(result.graph_data_name.c_str(), name_size);
    // write time
    ofs.write(reinterpret_cast<const char *>(&result.time), sizeof(result.time));
    ofs.close();
}

inline void GEDResultToBinary(const std::string &output_path, std::vector<GEDEvaluation> &results) {
    // set file path
    if (results.empty()) {
        std::cerr << "No results to save." << std::endl;
        return;
    }
    const std::string graph_data_name = results[0].graph_data_name;
    const std::string file_path = output_path + graph_data_name + "_ged_mapping.bin";
    // create output target_path if it does not exist
    if (!std::filesystem::exists(output_path)) {
        std::filesystem::create_directory(output_path);
    }
    // open binary file
    std::ofstream ofs(file_path, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for writing: " << file_path << std::endl;
        return;
    }
    for (const auto& result : results) {
        int source_id = result.graph_ids.first;
        int target_id = result.graph_ids.second;
        const std::vector<NodeId>& i_to_j = result.node_mapping.first;
        const std::vector<NodeId>& j_to_i = result.node_mapping.second;
        // write source_id and target_id
        ofs.write(reinterpret_cast<const char *>(&source_id), sizeof(source_id));
        ofs.write(reinterpret_cast<const char *>(&target_id), sizeof(target_id));
        // write size of i_to_j
        size_t size_i_to_j = i_to_j.size();
        ofs.write(reinterpret_cast<const char *>(&size_i_to_j), sizeof(size_i_to_j));
        // write i_to_j
        ofs.write(reinterpret_cast<const char *>(i_to_j.data()), size_i_to_j * sizeof(NodeId));
        // write size of j_to_i
        size_t size_j_to_i = j_to_i.size();
        ofs.write(reinterpret_cast<const char *>(&size_j_to_i), sizeof(size_j_to_i));
        // write j_to_i
        ofs.write(reinterpret_cast<const char *>(j_to_i.data()), size_j_to_i * sizeof(NodeId));
        // write graph_database_name
        size_t name_size = result.graph_data_name.size();
        ofs.write(reinterpret_cast<const char *>(&name_size), sizeof(name_size));
        ofs.write(result.graph_data_name.c_str(), name_size);
        // write time
        ofs.write(reinterpret_cast<const char *>(&result.time), sizeof(result.time));
    }
    ofs.close();

}

inline void BinaryToGEDResult(const std::string &input_path, GraphData<GraphStruct>& graph_data, GEDEvaluation &result) {
    // open binary file
    std::ifstream ifs(input_path, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error opening file for reading: " << input_path << std::endl;
        return;
    }
    // read source_id and target_id
    ifs.read(reinterpret_cast<char *>(&result.graph_ids.first), sizeof(result.graph_ids.first));
    ifs.read(reinterpret_cast<char *>(&result.graph_ids.second), sizeof(result.graph_ids.second));
    // read size of i_to_j
    size_t size_i_to_j;
    ifs.read(reinterpret_cast<char *>(&size_i_to_j), sizeof(size_i_to_j));
    // read i_to_j
    result.node_mapping.first.resize(size_i_to_j);
    ifs.read(reinterpret_cast<char *>(result.node_mapping.first.data()), size_i_to_j * sizeof(NodeId));
    // read size of j_to_i
    size_t size_j_to_i;
    ifs.read(reinterpret_cast<char *>(&size_j_to_i), sizeof(size_j_to_i));
    // read j_to_i
    result.node_mapping.second.resize(size_j_to_i);
    ifs.read(reinterpret_cast<char *>(result.node_mapping.second.data()), size_j_to_i * sizeof(NodeId));
    // read graph_database_name
    size_t name_size;
    ifs.read(reinterpret_cast<char *>(&name_size), sizeof(name_size));
    result.graph_data_name.resize(name_size);
    ifs.read(&result.graph_data_name[0], name_size);
    // read time
    ifs.read(reinterpret_cast<char *>(&result.time), sizeof(result.time));
    // set graphs
    result.graphs = {graph_data[result.graph_ids.first], graph_data[result.graph_ids.second]};
}

inline void BinaryToGEDResult(const std::string &input_path, GraphData<GraphStruct>& graph_data, std::vector<GEDEvaluation> &results) {
    // open binary file
    std::ifstream ifs(input_path, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error opening file for reading: " << input_path << std::endl;
        return;
    }
    while (ifs.peek() != EOF) {
        GEDEvaluation result;
        // read source_id and target_id
        ifs.read(reinterpret_cast<char *>(&result.graph_ids.first), sizeof(result.graph_ids.first));
        ifs.read(reinterpret_cast<char *>(&result.graph_ids.second), sizeof(result.graph_ids.second));
        // read size of i_to_j
        size_t size_i_to_j;
        ifs.read(reinterpret_cast<char *>(&size_i_to_j), sizeof(size_i_to_j));
        // read i_to_j
        result.node_mapping.first.resize(size_i_to_j);
        ifs.read(reinterpret_cast<char *>(result.node_mapping.first.data()), size_i_to_j * sizeof(NodeId));
        // read size of j_to_i
        size_t size_j_to_i;
        ifs.read(reinterpret_cast<char *>(&size_j_to_i), sizeof(size_j_to_i));
        // read j_to_i
        result.node_mapping.second.resize(size_j_to_i);
        ifs.read(reinterpret_cast<char *>(result.node_mapping.second.data()), size_j_to_i * sizeof(NodeId));
        // read graph_database_name
        size_t name_size;
        ifs.read(reinterpret_cast<char *>(&name_size), sizeof(name_size));
        result.graph_data_name.resize(name_size);
        ifs.read(&result.graph_data_name[0], name_size);
        // read time
        ifs.read(reinterpret_cast<char *>(&result.time), sizeof(result.time));
        // set graphs
        result.graphs = {graph_data[result.graph_ids.first], graph_data[result.graph_ids.second]};
        // push back result
        results.push_back(result);
    }
    ifs.close();
}

inline void MergeBinaries(const std::string& input_path, GraphData<GraphStruct>& graph_data, std::vector<GEDEvaluation> &results) {
    // read all mappings from binary files in the input path
    for (const auto& entry : std::filesystem::directory_iterator(input_path)) {
        if (entry.path().extension() == ".bin" && entry.path().filename().string().find("_ged_mapping") != std::string::npos) {
            results.emplace_back();
            BinaryToGEDResult(entry.path().string(), graph_data, results.back());
        }
    }
}




#endif //GEDEXAMPLE_GEDFUNCTIONS_H