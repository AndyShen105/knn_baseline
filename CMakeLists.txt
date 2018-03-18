cmake_minimum_required(VERSION 3.9)
project(Baseline)

set(CMAKE_CXX_STANDARD 11)

add_executable(Baseline main.cpp data.cpp main.h knn.cpp knn.h sort.h sort.cpp matrixMulCPU.h matrixMulCPU.cu)