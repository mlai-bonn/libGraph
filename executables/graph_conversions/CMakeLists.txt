cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 20)
project(graph_conversions)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(../../include)
link_directories(../../include)

add_executable(ConvertGraphs convert_graphs.cpp ../../include/libGraph.h)
