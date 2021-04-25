cmake_minimum_required(VERSION 3.17)
project(GraphCut)

set(CMAKE_CXX_STANDARD 17)

add_executable(graph_cut main.cpp stb/stb_lib.cpp)
add_executable(dft_test dft_test.cpp stb/stb_lib.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(graph_cut PUBLIC OpenMP::OpenMP_CXX)
endif()