cmake_minimum_required(VERSION 3.20)
project(Google_tests)
set(CMAKE_CXX_STANDARD 20)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(../include)
link_directories(../include)

# adding the Google_Tests_run target
add_executable(Google_Tests_run TestsMain.cpp ../include/libGraph.h)
target_link_libraries(Google_Tests_run gtest gtest_main)