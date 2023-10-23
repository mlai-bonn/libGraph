cmake_minimum_required(VERSION 3.20)
project(Google_tests)
set(CMAKE_CXX_STANDARD 20)
add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

include_directories(../include)
link_directories(../include)

# adding the Google_Tests_run target
add_executable(Google_Tests_run GraphAlgorithmsTest.cpp ../include/libGraph.h)

# linking Google_Tests_run with DateConverter_lib which will be tested
#target_link_libraries(Google_Tests_run include)

target_link_libraries(Google_Tests_run gtest gtest_main)