project(knn)

set(CMAKE_CXX_STANDARD 11)
#add_executable(knn data.cpp lsh.cpp lsh.h knn.cpp knn.h  v_lsh.cpp v_lsh.h main.cpp main.h util.cpp util.h hashFunc.cpp hashFunc.h)
#add_executable(testhash hashFunc.cpp hashFunc.h data.cpp)
#add_executable(v_lsh v_lsh.cpp data.cpp v_lsh.h)
#INCLUDE_DIRECTORIES(/home/andyshen/tools/gperftools/include)
#FIND_LIBRARY(PTHREAD_LIB NAMES pthread PATHS "/usr/lib/x86_64-linux-gnu/")
#MESSAGE(STATUS "Library path PTHREAD_LIB is " ${PTHREAD_LIB})
#FIND_LIBRARY(UNWIND_LIB NAMES unwind PATHS "/home/andyshen/tools/libunwind/lib/")
#MESSAGE(STATUS "Library path UNWIND_LIB is " ${UNWIND_LIB})
#FIND_LIBRARY(PROFILER_LIB NAMES profiler PATHS "/home/andyshen/tools/gperftools/lib/")
#MESSAGE(STATUS "Library path PROFILER_LIB is " ${PROFILER_LIB})
#
#TARGET_LINK_LIBRARIES(knn ${PTHREAD_LIB})
#TARGET_LINK_LIBRARIES(knn ${UNWIND_LIB})
#TARGET_LINK_LIBRARIES(knn ${PROFILER_LIB})
#
#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
cmake_minimum_required(VERSION 3.12)
set(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda-9.0)
SET( CMAKE_CXX_FLAGS "-std=c++11")
 #找CUDA
find_package(CUDA REQUIRED)

 #设置公共源文件
file(GLOB COMMON common.h common.cpp)

# 添加*.cu源文件
file(GLOB CU main.cu)

 #添加包含目录
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(/usr/local/cuda-9.0/include)
include_directories(/usr/local/cuda-9.0/samples/common/inc)


 #添加可执行文件
CUDA_ADD_EXECUTABLE(First # target名
        main.cu  # 添加*.cu源文件
        data.cpp
        data.h
        lsh.cpp
        lsh.h
        v_lsh.cpp
        v_lsh.h
        util.h
        hashFunc.cpp
        hashFunc.h
        )
