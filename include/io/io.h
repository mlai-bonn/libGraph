#include <vector>
#include <fstream>
#include <sstream>
#include <ostream>

#ifndef IO_IO_H
#define IO_IO_H
template <typename T>
class DataIO{
public:
    static void WriteTrivial(const std::string& s, const std::vector<T>& data );
    static void ReadTrivial(const std::string& s, std::vector<T>& data );
    static void WriteTrivialMatrix(const std::string& s, const std::vector<std::vector<T>>& data );
    static void ReadTrivialMatrix(const std::string& s, std::vector<std::vector<T>>& data );
    static void ToLatexTable(const std::string& path, const std::vector<std::vector<T>>& data);
};
#include "io.txx"

static void ArgumentRange(char** argv, int argc, int i, std::vector<int>& output){
    output.clear();
    int min = std::stoi(argv[i+1]);
    int max;
    int step = 1;
    if (i+2 < argc && std::string(argv[i+2]).find("..") != std::string::npos){
        if (std::string(argv[i+3]).find("--") == std::string::npos) {
            max = std::stoi(argv[i+3]);
            if (i + 4 < argc && std::string(argv[i + 4]).find("..") != std::string::npos && std::string(argv[i+5]).find("--") == std::string::npos) {
                step = std::stoi(argv[i+5]);
            }
            output.resize((max - min + step)/step);
            min -= step;
            std::generate(output.begin(), output.end(), [&min, step] {return min += step;});
        }
    }
    else {
        for (int j = i + 1; j < argc; ++j) {
            if (std::string(argv[j]).find("--") != std::string::npos) {
                break;
            } else {
                output.emplace_back(std::stoi(argv[j]));
            }
        }
    }
    std::cout << min << " " << max << " " << step << std::endl;
}
static void ArgumentList(char** argv, int argc, int i, std::vector<std::string>& output) {
    output.clear();
    for (int j = i + 1; j < argc; ++j) {
        if (std::string(argv[j]).find("--") != std::string::npos) {
            break;
        } else {
            output.emplace_back(argv[j]);
        }
    }
}

static void ArgumentList(char** argv, int argc, int i, std::vector<int>& output) {
    output.clear();
    for (int j = i + 1; j < argc; ++j) {
        if (std::string(argv[j]).find("--") != std::string::npos) {
            break;
        } else {
            output.emplace_back(std::stoi(argv[j]));
        }
    }
}


#endif
