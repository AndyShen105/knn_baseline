project(knn_baseline)

set(CMAKE_CXX_STANDARD 11)

#add_executable(knn main.cpp data.cpp main.h knn.cpp knn.h sort.h sort.cpp)
#add_executable(v_lsh v_lsh.cpp data.cpp v_lsh.h)
add_executable(knn  data.cpp lsh.cpp lsh.h knn.cpp knn.h util.h v_lsh.cpp v_lsh.h main.cpp main.h util.cpp util.h)