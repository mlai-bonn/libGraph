cmake_minimum_required(VERSION 3.23)
project(testGraphLib)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(Test test.cpp load_save.h graph_functions.h ../../include/libGraph.h ../../include/Hops/Hops.h ../../include/Hops/RootedPattern.h hops.h ../../include/io/FileEvaluation.h ../../include/io/io.h ../../include/io/io.txx ../../include/io/StaticFunctions.h)
