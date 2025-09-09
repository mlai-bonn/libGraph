cmake_minimum_required(VERSION 3.20)
project(GEDExample)
set(CMAKE_CXX_STANDARD 20)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(../../../include)
link_directories(../../../include)

# add ged library
# set gedlib root to ../../../gedlib
# Add the following directories to your include directories:
#
#    <GEDLIB_ROOT>
#    <GEDLIB_ROOT>/ext/boost.1.89.0
#    <GEDLIB_ROOT>/ext/eigen.3.3.4/Eigen
#    <GEDLIB_ROOT>/ext/nomad.3.8.1/src
#    <GEDLIB_ROOT>/ext/nomad.3.8.1/ext/sgtelib/src
#    <GEDLIB_ROOT>/ext/lsape.5/include
#    <GEDLIB_ROOT>/ext/libsvm.3.22
#    <GEDLIB_ROOT>/ext/fann.2.2.0/include

# Add the following directories to your link directories:
#
#    <GEDLIB_ROOT>/ext/nomad.3.8.1/lib
#    <GEDLIB_ROOT>/ext/libsvm.3.22
#    <GEDLIB_ROOT>/ext/fann.2.2.0/lib

set(GEDLIB_ROOT ../../../external/gedlib)
include_directories(${GEDLIB_ROOT} ${GEDLIB_ROOT}/ext/boost.1.89.0 ${GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen ${GEDLIB_ROOT}/ext/nomad.3.8.1/src ${GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src ${GEDLIB_ROOT}/ext/lsape.5/cpp/include ${GEDLIB_ROOT}/ext/libsvm.3.22 ${GEDLIB_ROOT}/ext/fann.2.2.0/include)
link_directories(${GEDLIB_ROOT}/ext/nomad.3.8.1/lib ${GEDLIB_ROOT}/ext/libsvm.3.22 ${GEDLIB_ROOT}/ext/fann.2.2.0/lib )


add_executable(GEDExample main.cpp ../../../include/libGraph.h ${GEDLIB_ROOT}/src/env/ged_env.hpp)

#Link your application against the following shared libraries:
#
#    <GEDLIB_ROOT>/ext/fann.2.2.0/lib/libdoublefann.2.dylib
#    <GEDLIB_ROOT>/ext/libsvm.3.22/libsvm.so
#    <GEDLIB_ROOT>/ext/nomad.3.8.1/libnomad.so

target_link_libraries(GEDExample libsvm.so libnomad.so libdoublefann.so.2)