cmake_minimum_required(VERSION 3.17)
project(GraphCut)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 17)

# library flag
set(CMAKE_CXX_FLAGS "-lfftw3")

# include directories
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

# executables
add_executable(graph_cut_baseline main.cpp stb/stb_lib.cpp)
add_executable(graph_cut_debug main.cpp stb/stb_lib.cpp)
add_executable(graph_cut main.cpp stb/stb_lib.cpp)

# target definitions
target_compile_definitions(graph_cut_baseline PUBLIC)
target_compile_definitions(graph_cut_debug PUBLIC DEBUG USE_FFTW)
target_compile_definitions(graph_cut PUBLIC USE_FFTW)

# link libraries
target_link_libraries(graph_cut_baseline fftw3)
target_link_libraries(graph_cut_debug fftw3)
target_link_libraries(graph_cut fftw3)
