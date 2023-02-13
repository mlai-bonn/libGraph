//
// Created by florian on 21.10.21.
//

#ifndef CLOSURES_FILEEVALUATION_H
#define CLOSURES_FILEEVALUATION_H


#include <string>
#include <vector>
#include <iomanip>
#include <unordered_map>

class FileEvaluation {
    std::string out_path;
    std::string name;
    std::vector<std::string> headers_all;
    std::unordered_map<std::string, std::vector<std::string>> values_all;
    std::vector<std::string> headers_summary;
    std::unordered_map<std::string, std::vector<std::string>> values_summary;
    std::string extension = ".csv";
public:
    FileEvaluation()= default;
    explicit FileEvaluation(const std::string& out_path, const std::string& eval_name = "default", const std::string& extension = ".csv") : out_path(out_path), name(eval_name), extension(extension){};
    void save(bool summary = false, bool both = true, std::_Ios_Openmode mode = std::ios_base::app);
    void headerValueInsert(const std::vector<std::string>& new_header, const std::vector<std::string>& new_values, int insert_position=-1, bool summary = false, bool both=false);
    void clear();

};


inline void FileEvaluation::save(bool summary, bool both, std::_Ios_Openmode mode) {
    if (!out_path.empty()) {
        if (summary && both) {
            save(false, false, mode);
        }
        std::string appendix = extension;
        if (summary) {
            appendix = "_summary" + extension;
        }
        bool newFile = !std::filesystem::exists(this->out_path + this->name + appendix);
        std::ofstream fs;
        fs.open(this->out_path + this->name + appendix, mode);

        // write the file headers
        std::vector<std::string> *headers = &headers_all;
        std::unordered_map<std::string, std::vector<std::string>> *values = &values_all;

        if (summary) {
            headers = &headers_summary;
            values = &values_summary;
        }
        if (newFile) {
            for (const auto &header: *headers) {
                fs << "," << header;
            }
            fs << std::endl;
        }


        int max_size = 0;
        for (const auto &header: *headers) {
            max_size = std::max(max_size, (int) (*values)[header].size());
        }

        fs << std::fixed;
        for (int i = 0; i < max_size; ++i) {
            for (auto const& header : *headers) {
                if ((*values)[header].size() <= i){
                    fs << "," << (*values)[header][(*values)[header].size() - 1];
                }
                else {
                    fs << "," << (*values)[header][i];
                }
            }
            fs << std::endl;
        }
        fs << std::scientific;
        fs.close();
    }
}

inline void FileEvaluation::headerValueInsert(const std::vector<std::string> &new_header, const std::vector<std::string> &new_values, int insert_position, bool summary, bool both) {
    if (new_header.size() != new_values.size()){
        throw std::length_error("Header and Value lengths do not fit!");
    }

    for (int i=0; i<new_header.size(); ++i) {
        const std::string header = new_header[i];

        if (!summary || both) {
            if (std::find(headers_all.begin(), headers_all.end(), header) == headers_all.end()) {
                headers_all.emplace_back(header);
                values_all.emplace(
                        std::pair<std::string, std::vector<std::string>>{header, std::vector<std::string>()});
            }

            std::vector<std::string> &vals = values_all[header];
            if (insert_position == -1){
                insert_position = (int) vals.size();
            }
            if (insert_position == -2) {
                vals.resize(1, new_values[i]);
            } else if (vals.size() <= insert_position) {
                vals.resize(insert_position + 1, "");
                vals[insert_position] = new_values[i];
            } else{
                vals[insert_position] = new_values[i];
            }
        }

        if (summary) {
            if (std::find(headers_summary.begin(), headers_summary.end(), header) == headers_summary.end()) {
                headers_summary.emplace_back(header);
                values_summary.emplace(
                        std::pair<std::string, std::vector<std::string>>{header, std::vector<std::string>()});
            }
            std::vector<std::string> &vals = values_summary[header];
            if (insert_position == -1){
                insert_position = (int) vals.size();
            }
            if (insert_position == -2) {
                vals.resize(1, new_values[i]);
            } else if (vals.size() <= insert_position) {
                vals.resize(insert_position + 1, "");
                vals[insert_position] = new_values[i];
            } else {
                vals[insert_position] = new_values[i];
            }
        }
    }
}

inline void FileEvaluation::clear() {
    headers_all.clear();
    headers_summary.clear();
    values_all.clear();
    values_summary.clear();
}



#endif //CLOSURES_FILEEVALUATION_H
