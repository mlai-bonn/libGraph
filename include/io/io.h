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

};
#include "io.txx"

#endif
