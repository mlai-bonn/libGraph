cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 20)
project(libGraph)
add_subdirectory(Google_tests)

add_executable(ConvertGraphs executables/convert_graphs.cpp /include/libGraph.h)


