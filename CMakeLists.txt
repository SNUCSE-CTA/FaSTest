cmake_minimum_required(VERSION 3.21)
project(Fastest)

message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -Ofast")
FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(GSL REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})
include_directories(lib)


add_executable(Fastest driver/subgraph-cardinality-estimation.cc)