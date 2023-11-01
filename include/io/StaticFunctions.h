//
// Created by Florian on 13.04.2021.
//

#ifndef IOLIB_STATICFUNCTIONS_H
#define IOLIB_STATICFUNCTIONS_H


#include <string>
#include <map>
#include <math.h>
#include <numeric>
#include <algorithm>
#include "iostream"
#include <random>
#include <filesystem>

class GraphParams;


class StaticFunctionsLib {
public:
    template<class T>
    static T mean(const std::vector<T>& vector);

    template<class T>
    static T standard_deviation(const std::vector<T>& vector);

    template<class T>
    static T median(std::vector<T>& vector);

    template<class T1, class T2>
    static std::string print(T1 &Object);

    template<class T1, class T2>
    static std::string print(T1 &Object, bool pair);

    template<class T2>
    static std::string pairToString(const T2 &object);

    // get the paths of directories inside path using the names in names
    static std::vector<std::string> directory_paths(const std::string & path, const std::vector<std::string>& names);
    // get the paths of files inside path using the names in names
    static std::vector<std::string> file_paths(const std::string & path, const std::vector<std::string>& names = {}, const std::vector<std::string>& extensions = {});

    static std::string printMap(const std::map<int, int> &map);
    static void write_csv(const std::string& path, std::vector<std::vector<std::string>>& data, const char& delimiter = ',');
    static void load_csv(const std::string &path, std::vector<std::vector<std::string>>& out, const char& delimiter = ',');
    static void load_csv(const std::string &path, std::map<std::string, std::vector<std::string>>& out, const char& delimiter = ',');

    static void get_k_from_n(std::vector<int>& output, int k, int n, int seed = 0, std::vector<int>* input = nullptr);

    static size_t get_vector_position(const std::vector<size_t>& vector_lengths, const std::vector<int>& positions);

    // Comparison function to sort the vector elements
    // by second element of tuples
    static bool sortbysecond(const std::tuple<std::string, int>& a,
                      const std::tuple<std::string, int>& b)
    {
        return (std::get<1>(a) < std::get<1>(b));
    }

    static auto above(double threshold) {
        // This captures a copy of threshold
        return [=](double value) {
            return value >= threshold;
        };
    };

    static void save(const std::string &path, const std::vector<double> &values, const std::string &extension);

    static void load(const std::string &path, std::vector<double> &values);

    static void PrintStream(std::stringstream& stringstream);

    static void GetkFoldIndices(int kFold, int dataSize, std::vector<std::pair<std::vector<int>, std::vector<int>>>& indices, int seed);

    template<class T1, class T2>
    static std::string mapToString(std::unordered_map<T1, T2> &inputMap);

    template<class T1, class T2>
    static std::string mapToString(std::map<T1, T2> &inputMap);

    template<class T1, class T2>
    static void print(std::unordered_map<T1, T2> &inputMap);

    template<class T1, class T2>
    static void print(std::map<T1, T2> &inputMap);

    template<class T1>
    static std::string vectorToString(std::vector<T1> &inputVector);
    template<class T1>
    static void stringToVector(std::vector<T1>& out, std::string& string);

    template<class T1>
    static void print(std::vector<T1> &inputVector);

    static void saveValuesToFile(const std::string &path, const std::vector<std::string> &header,
                          const std::vector<std::string> &values, std::_Ios_Openmode mode);
};

template <class T1, class T2>
inline std::string StaticFunctionsLib::print(T1& Object) {


    return "{" + std::accumulate(std::begin(Object),
                           std::end(Object),
                           std::string{},
                           [](const std::string &a, const T2 &b) {
                               return a.empty() ? std::to_string(b)
                                                : a + ", " + std::to_string(b);
                           }) + "}";
}

template<class T1, class T2>
inline std::string StaticFunctionsLib::print(T1 &Object, bool pair) {
    if (pair) {
        return std::accumulate(std::begin(Object),
                               std::end(Object),
                               std::string{},
                               [](const std::string &a, const T2 &b) {
                                   return a.empty() ? '"' + pairToString(b) + '"'
                                                    : a + ", " + '"' + pairToString(b) + '"';
                               });
    }
    return "";
}

template<class T2>
inline std::string StaticFunctionsLib::pairToString(const T2 &object) {
    return "(" + std::to_string(object.first) + "," + std::to_string(object.second) + ")";
}

template<class T>
inline T StaticFunctionsLib::mean(const std::vector<T> &vector) {
    return std::accumulate(vector.begin(),  vector.end(), 0.0)/(double) vector.size();
}

template<class T>
inline T StaticFunctionsLib::standard_deviation(const std::vector<T> &vector) {
    T m = StaticFunctionsLib::mean(vector);
    double accum = 0.0;
    std::for_each (vector.begin(),  vector.end(), [&](const double d) {
        accum += (d - m) * (d - m);
    });
    return std::sqrt(accum / (vector.size()-1));
}

template<class T>
inline T StaticFunctionsLib::median(std::vector<T> &vector) {
    if (vector.size() % 2 == 0) {
        const auto median_it1 = vector.begin() + vector.size() / 2 - 1;
        const auto median_it2 = vector.begin() + vector.size() / 2;

        std::nth_element(vector.begin(), median_it1 , vector.end());
        const auto e1 = *median_it1;

        //std::nth_element(vector.begin(), median_it2 , vector.end());
        const auto e2 = *median_it2;

        return (e1 + e2) / 2;

    } else {
        const auto median_it = vector.begin() + vector.size() / 2;
        //std::nth_element(vector.begin(), median_it , vector.end());
        return *median_it;
    }
}

inline std::string StaticFunctionsLib::printMap(const std::map<int, int> &map) {
    std::string out = "{";
    for (auto const & [key, value] : map) {
        out += "(";
        out += std::to_string(key);
        out += ":";
        out += std::to_string(value);
        out += ")";
    }
    out += "}";
    return out;
}



inline void StaticFunctionsLib::save(const std::string &path,const std::vector<double> &values, const std::string& extension) {
    std::string complete_path = path + extension;

    if (extension.find("bin") != std::string::npos){
        std::ofstream file(complete_path, std::ios::out | std::ios::binary);
        size_t size = (values.size());
        file.write((char*) (&size), sizeof(size_t));
        for (auto const val: values) {
            file.write((char*) (&val), sizeof(double));
        }
        file.close();
    }
    else {
        std::ofstream file(complete_path, std::ios::out);
        for (auto const val: values) {
            file << val << std::endl;
        }
        file.close();
    }
}
inline void StaticFunctionsLib::load(const std::string &path, std::vector<double> &values) {
    std::string extension = std::filesystem::path(path).extension();
    if (extension.find("bin") != std::string::npos){
        std::ifstream file(path, std::ios::in | std::ios::binary);
        size_t size = 0;
        file.read((char*) (&size), sizeof(size_t));
        values = std::vector<double>(size, 0);
        for (auto & val: values) {
            file.read((char*) (&val), sizeof(double));
        }
        file.close();
    }
    else {
        std::ifstream file(path, std::ios::in);
        std::string row, item;
        values.clear();
        while(std::getline(file, row)) {
            std::istringstream ss(row);
            while (std::getline(ss, item)) {
                values.emplace_back(std::stod(item));
            }
        }
        file.close();
    }
}

inline void StaticFunctionsLib::write_csv(const std::string& path, std::vector<std::vector<std::string>>& data, const char& delimiter)
{
    std::ofstream file;
    file.open(path);
    for (auto const& row : data) {
        for (int i = 0; i < row.size(); ++i) {
            if (i != row.size() - 1){
                file << row[i] << ",";
            }
            else{
                file << row[i];
            }
        }
        file << std::endl;
    }
    file.close();
}

inline void StaticFunctionsLib::load_csv(const std::string &path, std::vector<std::vector<std::string>>& out, const char& delimiter) {
    std::string row, item;
    std::ifstream in(path);
    std::vector<std::string> R;
    while(std::getline(in, row))
    {
        R.clear();
        std::istringstream ss(row);
        while (std::getline(ss, item, delimiter)) {
            R.push_back(item);
        }
        out.push_back( R );
    }
    in.close();
}

inline void StaticFunctionsLib::load_csv(const std::string &path, std::map<std::string, std::vector<std::string>>& out, const char& delimiter) {
    std::string row, item;
    std::ifstream in(path);
    std::vector<std::string> R;
    std::vector<std::string> col_names;
    int row_counter = 0;
    while(std::getline(in, row))
    {
        R.clear();
        std::istringstream ss(row);
        int col_counter = 0;
        while (std::getline(ss, item, delimiter)) {
            if (row_counter == 0) {
                col_names.push_back(item);
                out[item] = std::vector<std::string>();
            }
            else{
                out[col_names[col_counter]].push_back(item);
            }
            ++col_counter;
        }
        ++row_counter;
    }
    in.close();
}

inline void StaticFunctionsLib::PrintStream(std::stringstream& stringstream) {
    std::cout << stringstream.str();
}

inline void StaticFunctionsLib::get_k_from_n(std::vector<int> &output, int k, int n, int seed, std::vector<int> *input) {
    std::mt19937_64 gen(seed);
    if (input != nullptr && input->size() == n){
        for (int i = 0; i < k; ++i) {
            int rand = std::uniform_int_distribution<int>(i, n-1)(gen);
            output.emplace_back((*input)[rand]);
            std::swap((*input)[rand], (*input)[i]);
        }
    }
    else{
        std::vector<int> vec = std::vector<int>(n, 0);
        std::iota(vec.begin(), vec.end(), 0);
        output.clear();
        for (int i = 0; i < k; ++i) {
            int rand = std::uniform_int_distribution<int>(i, n-1)(gen);
            output.emplace_back(vec[rand]);
            std::swap(vec[rand], vec[i]);
        }
    }
}

inline size_t StaticFunctionsLib::get_vector_position(const std::vector<size_t>& vector_lengths, const std::vector<int>& positions) {
    if (vector_lengths.size() != positions.size()){
        throw std::range_error("The lengths of the two vectors have to be equal.");
    }
    size_t position = 0;
    for (int i = 0; i < vector_lengths.size(); ++i) {
        size_t result = 1;
        for (int j = i+1; j < vector_lengths.size(); ++j) {
            result *= vector_lengths[j];
        }
        position += positions[i] * result;
    }
    return position;
}

inline void StaticFunctionsLib::GetkFoldIndices(int kFold, int dataSize, std::vector<std::pair<std::vector<int>, std::vector<int>>>& indices, int seed) {
    std::vector<int> ind = std::vector<int>(dataSize);
    std::iota(ind.begin(), ind.end(), 0);
    std::mt19937_64 gen(0);
    int maxIndex = (int) ind.size() - 1;
    for (int k = 0; k < kFold; ++k) {
        indices.emplace_back();
        int size = (int) std::round((float) dataSize/(float) kFold);
        int counter = 0;
        while(maxIndex >= 0 && counter < size){
            int rand = std::uniform_int_distribution<int>(0, maxIndex)(gen);
            indices[k].first.emplace_back(ind[rand]);
            std::swap(ind[rand], ind[maxIndex]);
            --maxIndex;
            ++counter;
        }
    }
    for (int k = 0; k < kFold; ++k) {
        for (auto x : ind) {
            if (std::find(indices[k].first.begin(), indices[k].first.end(), x) == indices[k].first.end()){
                indices[k].second.emplace_back(x);
            }
        }
    }
}

template<typename T1, typename T2>
inline std::string StaticFunctionsLib::mapToString(std::unordered_map<T1, T2> &inputMap) {
    std::stringstream ss;
    ss << std::fixed << "{";
    for (const auto &[k, v] : inputMap) {
        if (v == (*inputMap.begin()).second){
            ss << k << " : " << v;
        }
        else {
            ss << " " << k << " : " << v;
        }
    }
    ss << "}";
    return ss.str();
}
template<typename T1, typename T2>
inline std::string StaticFunctionsLib::mapToString(std::map<T1, T2> &inputMap) {
    std::stringstream ss;
    ss << std::fixed << "{";
    for (const auto &[k, v] : inputMap) {
        if (v == (*inputMap.begin()).second){
            ss << k << " : " << v;
        }
        else {
            ss << " " << k << " : " << v;
        }
    }
    ss << "}";
    return ss.str();
}
template<typename T1, typename T2>
inline void StaticFunctionsLib::print(std::unordered_map<T1, T2> &inputMap) {
    std::cout << mapToString<T1, T2>(inputMap) << std::endl;
}
template<typename T1, typename T2>
inline void StaticFunctionsLib::print(std::map<T1, T2> &inputMap) {
    std::cout << mapToString<T1, T2>(inputMap) << std::endl;
}
template<typename T1>
inline std::string StaticFunctionsLib::vectorToString(std::vector<T1> &inputVector) {
    std::stringstream ss;
    ss << std::fixed << "[";
    for (const auto &v : inputVector) {
        if (v == *inputVector.begin()){
            ss << v;
        }
        else{
            ss << " " << v;
        }
    }
    ss << "]";
    return ss.str();
}

template<typename T1>
inline void StaticFunctionsLib::print(std::vector<T1> &inputVector) {
    std::cout << vectorToString<T1>(inputVector) << std::endl;
}

inline void StaticFunctionsLib::saveValuesToFile(const std::string& path, const std::vector<std::string>& header, const std::vector<std::string>& values, std::_Ios_Openmode mode) {
    bool newFile = std::filesystem::exists(path);
    std::ofstream fs;
    fs.open(path, mode);
    // write the file headers

    if (!newFile) {
        for (const auto& string : header) {
            fs << "," << string;
        }
        fs << std::endl;
    }
    fs << std::fixed;
    for (const auto& string : values) {
        fs << "," << string;
    }
    fs << std::endl;
    fs << std::scientific;
    fs.close();
}

inline std::vector<std::string> StaticFunctionsLib::directory_paths(const std::string& path, const std::vector<std::string>& names) {
    // search path for directories that match with entries in names
    std::vector<std::string> paths;
    // iterate recursively over the directory
    for (const auto& entry : std::filesystem::recursive_directory_iterator(path)) {
        // check if the entry is a directory
        if (std::filesystem::is_directory(entry)) {
            // check if the directory name is in the list of names
            if (std::find(names.begin(), names.end(), entry.path().filename()) != names.end()) {
                paths.emplace_back(entry.path());
            }
        }
    }
    // return the paths
    return paths;
}

inline std::vector<std::string>
StaticFunctionsLib::file_paths(const std::string &path, const std::vector<std::string> &names, const std::vector<std::string> &extensions) {
    // search path for files that match with entries in names
    std::vector<std::string> paths;
    // iterate recursively over the directory
    for (const auto& entry : std::filesystem::recursive_directory_iterator(path)) {
        // get entry
        std::string entry_path = entry.path();


                // check if the entry is a regular file
                if (std::filesystem::is_regular_file(entry)) {
                    // get filename
                    std::string filename = entry.path().filename();
                    // get file extension
                    std::string extension = entry.path().extension();
                    // check if the file extensions is in the list of extensions
                    if (!extensions.empty()) {
                        if (std::find(extensions.begin(), extensions.end(), entry.path().extension()) !=
                            extensions.end()) {
                            // check if the file name is in the list of names
                            paths.emplace_back(entry.path());
                        }
                    }
                    // check if the file name is in the list of names
                    if (!names.empty()) {
                        if (std::find(names.begin(), names.end(), entry.path().filename()) != names.end()) {
                            paths.emplace_back(entry.path());
                        }
                    }
                    if (names.empty() && extensions.empty()) {
                        paths.emplace_back(entry.path());
                    }
            }
        }
    // return the paths
    return paths;
}


template<class T1>
inline void StaticFunctionsLib::stringToVector(std::vector<T1>& out, std::string &string) {
    string.erase(remove(string.begin(), string.end(), '['), string.end());
    string.erase(remove(string.begin(), string.end(), ']'), string.end());

    std::stringstream ss(string);
    std::string segment;
    std::vector<std::string> list;

    while(std::getline(ss, segment, ' '))
    {
        list.emplace_back(segment);
    }
    for (const auto& x : list) {
        out.emplace_back((T1)x);
    }
}

#endif //IOLIB_STATICFUNCTIONS_H
