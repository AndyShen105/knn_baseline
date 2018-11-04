project(knn)

set(CMAKE_CXX_STANDARD 11)
add_executable(knn data.cpp lsh.cpp lsh.h knn.cpp knn.h util.h v_lsh.cpp v_lsh.h main.cpp main.h util.cpp util.h)
#add_executable(knn main.cpp data.cpp main.h knn.cpp knn.h sort.h sort.cpp)
#add_executable(v_lsh v_lsh.cpp data.cpp v_lsh.h)
INCLUDE_DIRECTORIES(/home/andyshen/tools/gperftools/include)
FIND_LIBRARY(PTHREAD_LIB NAMES pthread PATHS "/usr/lib/x86_64-linux-gnu/")
MESSAGE(STATUS "Library path PTHREAD_LIB is " ${PTHREAD_LIB})
FIND_LIBRARY(UNWIND_LIB NAMES unwind PATHS "/home/andyshen/tools/libunwind/lib/")
MESSAGE(STATUS "Library path UNWIND_LIB is " ${UNWIND_LIB})
FIND_LIBRARY(PROFILER_LIB NAMES profiler PATHS "/home/andyshen/tools/gperftools/lib/")
MESSAGE(STATUS "Library path PROFILER_LIB is " ${PROFILER_LIB})

TARGET_LINK_LIBRARIES(knn ${PTHREAD_LIB})
TARGET_LINK_LIBRARIES(knn ${UNWIND_LIB})
TARGET_LINK_LIBRARIES(knn ${PROFILER_LIB})

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

