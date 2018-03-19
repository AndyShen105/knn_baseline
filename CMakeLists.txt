cmake_minimum_required(VERSION 3.9)
project(Baseline)
include_directories(
        /Developer/NVIDIA/CUDA-8.0/include
        /Developer/NVIDIA/CUDA-8.0/samples/common/inc

)
set(CMAKE_CXX_STANDARD 11)
find_package(CUDA)
file(GLOB_RECURSE CURRENT_HEADERS  *.h *.hpp *.cuh)
file(GLOB CURRENT_SOURCES  *.cpp *.cu)
CUDA_ADD_EXECUTABLE(Baseline main.cpp data.cpp main.h knn.cpp knn.h sort.h sort.cpp matrixMulGPU.h matrixMulGPU.cu)