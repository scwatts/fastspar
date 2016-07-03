cmake_minimum_required(VERSION 3.5)
project(sparcpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -larmadillo")

set(SOURCE_FILES sparcpp.cpp)
add_executable(sparcpp ${SOURCE_FILES})
