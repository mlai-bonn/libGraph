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
template<typename T>
inline std::vector<T> CreateEditPath(const GEDEvaluation<T>& result)
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

    EditPath<T> edit_path;
    result.get_edit_path(edit_path, 0);
    // print edit path length
    std::cout << "Edit Path Length: " << edit_path.edit_path_graphs.size() - 1 << std::endl;
    std::vector<T> edit_path_graphs;
    for (const auto& g : edit_path.edit_path_graphs) {
        edit_path_graphs.emplace_back(g);
    }
    edit_path_graphs.back().set_name(edit_path.target_graph.GetName());
    // print all the graphs
    //for (const auto& g : edit_path_graphs.graphData) {
    //    std::cout << g << std::endl;
    //}
    return edit_path_graphs;
}

template<typename T>
inline void WriteEditPathInfo(const std::vector<GEDEvaluation<T>>& results, const GraphData<T>& graph_data, const std::string& edit_path_info) {
    // create also a MUTAG_edit_paths.bin storing source_id, step_id, target_id for each graph
    std::ofstream ofs(edit_path_info, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for writing: " << edit_path_info << std::endl;
        return;
    }
    for (const auto& result : results) {
        INDEX source_id = result.graph_ids.first;
        INDEX target_id = result.graph_ids.second;
        auto edit_path_graphs = CreateEditPath<DDataGraph>(result);
        for (INDEX step_id = 0; step_id < edit_path_graphs.size(); ++step_id) {
            // write source_id, step_id, target_id
            ofs.write(reinterpret_cast<const char *>(&source_id), sizeof(source_id));
            ofs.write(reinterpret_cast<const char *>(&step_id), sizeof(step_id));
            ofs.write(reinterpret_cast<const char *>(&target_id), sizeof(target_id));
        }
    }
    ofs.close();
}

inline void ReadEditPathInfo(std::string& edit_path_info, std::vector<std::tuple<INDEX, INDEX, INDEX>>& info) {
    // open binary file
    std::ifstream ifs(edit_path_info, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error opening file for reading: " << edit_path_info << std::endl;
        return;
    }
    info.clear();
    while (!ifs.eof()) {
        INDEX source_id;
        INDEX step_id;
        INDEX target_id;
        ifs.read(reinterpret_cast<char *>(&source_id), sizeof(source_id));
        ifs.read(reinterpret_cast<char *>(&step_id), sizeof(step_id));
        ifs.read(reinterpret_cast<char *>(&target_id), sizeof(target_id));
        info.emplace_back(std::make_tuple(source_id, step_id, target_id));
    }

}

template<typename T>
void CreateAllEditPaths(const std::vector<GEDEvaluation<T>> &results, const GraphData<T> &graph_data, const std::string &edit_path_output = "../Data/EditPaths/") {
    // check whether file already exists
    if (std::filesystem::exists(edit_path_output + graph_data.GetName() + "_edit_paths.bgf")) {
        std::cout << "Edit paths for " << graph_data.GetName() << " already exist." << std::endl;
        return;
    }
    GraphData<T> all_path_graphs;
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
        auto edit_path_graphs = CreateEditPath<T>(result);
        int path_counter = 0;
        for (auto& g : edit_path_graphs) {
            g.set_name(graph_data.GetName() + "_" + std::to_string(result.graph_ids.first) + "_" + std::to_string(result.graph_ids.second) + "_" + std::to_string(path_counter));
            all_path_graphs.add(g);
            ++path_counter;
        }
    }
    // save the final result in the tmp folder under datasetname_i.bgfs
    SaveParams params = {
        edit_path_output,
        graph_data.GetName() + "_edit_paths",
        GraphFormat::BGF,
        true,
    };
    all_path_graphs.Save(params);
    WriteEditPathInfo(results, graph_data, edit_path_output + graph_data.GetName() + "_edit_paths_info.bin");

}


template <typename T>
inline void GEDResultToBinary(const std::string &target_path, GEDEvaluation<T> &result) {
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
    // write distance
    ofs.write(reinterpret_cast<const char *>(&result.distance), sizeof(result.distance));
    // write lower bound
    ofs.write(reinterpret_cast<const char *>(&result.lower_bound), sizeof(result.lower_bound));
    // write upper bound
    ofs.write(reinterpret_cast<const char *>(&result.upper_bound), sizeof(result.upper_bound));
    ofs.close();
}

template <typename T>
inline void GEDResultToBinary(const std::string &output_path, std::vector<GEDEvaluation<T>> &results) {
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
        // write distance
        ofs.write(reinterpret_cast<const char *>(&result.distance), sizeof(result.distance));
        // write lower bound
        ofs.write(reinterpret_cast<const char *>(&result.lower_bound), sizeof(result.lower_bound));
        // write upper bound
        ofs.write(reinterpret_cast<const char *>(&result.upper_bound), sizeof(result.upper_bound));
    }
    ofs.close();
}


template<typename T>
inline void BinaryToGEDResult(const std::string &input_path, const GraphData<T>& graph_data, GEDEvaluation<T> &result) {
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
    // read distance
    ifs.read(reinterpret_cast<char *>(&result.distance), sizeof(result.distance));
    // read lower bound
    ifs.read(reinterpret_cast<char *>(&result.lower_bound), sizeof(result.lower_bound));
    // read upper bound
    ifs.read(reinterpret_cast<char *>(&result.upper_bound), sizeof(result.upper_bound));
    // set graphs
    result.graphs = {graph_data.graphData[result.graph_ids.first], graph_data.graphData[result.graph_ids.second]};
}

template<typename T>
void BinaryToGEDResult(const std::string &input_path, const GraphData<T>& graph_data, std::vector<GEDEvaluation<T>> &results) {
    // open binary file
    std::ifstream ifs(input_path, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error opening file for reading: " << input_path << std::endl;
        return;
    }
    while (ifs.peek() != EOF) {
        GEDEvaluation<T> result;
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
        // read distance
        ifs.read(reinterpret_cast<char *>(&result.distance), sizeof(result.distance));
        // read lower bound
        ifs.read(reinterpret_cast<char *>(&result.lower_bound), sizeof(result.lower_bound));
        // read upper bound
        ifs.read(reinterpret_cast<char *>(&result.upper_bound), sizeof(result.upper_bound));
        // set graphs
        result.graphs = {graph_data.graphData[result.graph_ids.first], graph_data.graphData[result.graph_ids.second]};
        // push back result
        results.push_back(result);
    }
    ifs.close();
}

template<typename T>
void MergeBinaries(const std::string& input_path, const std::string& search_string, const GraphData<T>& graph_data, std::vector<GEDEvaluation<T>> &results) {
    // read all mappings from binary files in the input path
    for (const auto& entry : std::filesystem::directory_iterator(input_path)) {
        if (entry.path().extension() == ".bin" && entry.path().filename().string().find(search_string) != std::string::npos) {
            results.emplace_back();
            BinaryToGEDResult(entry.path().string(), graph_data, results.back());
        }
    }
    // sort results by graph ids (first, second)
    std::ranges::sort(results, [](const GEDEvaluation<T>& a, const GEDEvaluation<T>& b) {
        if (a.graph_ids.first != b.graph_ids.first) {
            return a.graph_ids.first < b.graph_ids.first;
        }
        return a.graph_ids.second < b.graph_ids.second;
    });
}

inline void CSVFromGEDResults(const std::string &results_path, const std::vector<GEDEvaluation<DDataGraph>> &results) {
    // write distance csv file source_id, target_id, distance
    std::ofstream distance_ofs(results_path, std::ios::out);
    if (!distance_ofs) {
        std::cerr << "Error opening file for writing: " << results_path << std::endl;
        return;
    }
    // Write Header source_id, target_id, lower_bound, upper_bound, distance
    distance_ofs << "source_id, target_id, lower_bound, upper_bound, distance" << std::endl;
    for (const auto& result : results) {
        distance_ofs << result.graph_ids.first << "," << result.graph_ids.second << "," << result.lower_bound << "," << result.upper_bound << "," << result.distance << std::endl;
    }
    distance_ofs.close();
}

template<typename T>
void MergeGEDResults(const std::string &results_path, const std::string& search_string, const GraphData<T>& graph_data) {
    // Merge all mappings in tmp folder
    std::vector<GEDEvaluation<T>> merged_results;
    MergeBinaries(results_path + "tmp/", search_string, graph_data, merged_results);
    // Save final result
    GEDResultToBinary(results_path, merged_results);
    // remove all databasename related files in tmp folder
    for (const auto& entry : std::filesystem::directory_iterator(results_path + "tmp/")) {
        if (entry.path().string().find(graph_data.GetName() + "_") != std::string::npos) {
            std::filesystem::remove(entry.path());
        }
    }
}



#endif //GEDEXAMPLE_GEDFUNCTIONS_H