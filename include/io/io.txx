#include "io.h"
#include <ostream>
#include <fstream>
#include <iostream>

// This writes a vector of trivial data types.
template<typename T>
void DataIO<T>::WriteTrivial(const std::string& s, const std::vector<T>& data )
{
    std::ofstream out(s, std::ios::out | std::ios::binary);
    unsigned int len = data.size();
    out.write( (char*)&len, sizeof(len) );
    out.write( (const char*)&data[0], len * sizeof(T) );
    out.close();
}

// This reads a vector of trivial data types.
template<typename T>
void DataIO<T>::ReadTrivial(const std::string& s, std::vector<T>& data )
{
    std::ifstream in(s, std::ios::in | std::ios::binary);
    unsigned int len = 0;
    in.read( (char*)&len, sizeof(len) );
    data.resize(len);
    if( len > 0 ) in.read( (char*)&data[0], len * sizeof(T) );
    in.close();
}

template<typename T>
void DataIO<T>::WriteTrivialMatrix(const std::string& s, const std::vector<std::vector<T>>& data ){
    std::ofstream out(s, std::ios::out | std::ios::binary);
    unsigned int len = data.size();
    out.write( (char*)&len, sizeof(len) );
    for (auto const & dat : data) {
        len = dat.size();
        out.write( (char*)&len, sizeof(len) );
        out.write( (const char*)&dat[0], len * sizeof(T) );
    }
    out.close();
}
template<typename T>
void DataIO<T>::ReadTrivialMatrix(const std::string& s, std::vector<std::vector<T>>& data ){
    std::ifstream in(s, std::ios::in | std::ios::binary);
    unsigned int len = 0;
    unsigned int dat_len = 0;
    in.read( (char*)&len, sizeof(len) );
    data.resize(len);
    for (int i = 0; i < len; ++i) {
        auto dat = std::vector<T>();
        in.read( (char*)&dat_len, sizeof(dat_len) );
        dat.resize(dat_len);
        if( dat_len > 0 ) {
            in.read((char *) &dat[0], dat_len * sizeof(T));
            data[i] = dat;
        }
    }
    in.close();
}
