cmake_minimum_required(VERSION 3.17)
project(GraphCut)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 17)

# library flag
set(CMAKE_CXX_FLAGS "-lfftw3")

# include directories
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(graph_cut main.cpp stb/stb_lib.cpp)

# this for debug
add_executable(dft_test dft_test.cpp stb/stb_lib.cpp)

# link libraries
target_link_libraries(graph_cut fftw3)
target_link_libraries(dft_test fftw3)
