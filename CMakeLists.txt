cmake_minimum_required(VERSION 3.17)
project(GraphCut)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 17)

# library flag
set(CMAKE_CXX_FLAGS "-fopenmp -lfftw3_omp -lfftw3 -lm")

# include directories
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

# executables
add_executable(graph_cut_baseline main.cpp stb/stb_lib.cpp)
add_executable(graph_cut_fftw main.cpp stb/stb_lib.cpp)
add_executable(graph_cut_omp main.cpp stb/stb_lib.cpp)
add_executable(graph_cut_debug main.cpp stb/stb_lib.cpp)

# target definitions
target_compile_definitions(graph_cut_baseline PUBLIC)
target_compile_definitions(graph_cut_fftw PUBLIC USE_FFTW)
target_compile_definitions(graph_cut_omp PUBLIC USE_FFTW USE_OMP)
target_compile_definitions(graph_cut_debug PUBLIC DEBUG USE_FFTW USE_OMP)

# link libraries
target_link_libraries(graph_cut_baseline fftw3)
target_link_libraries(graph_cut_fftw fftw3)
target_link_libraries(graph_cut_omp fftw3)
target_link_libraries(graph_cut_debug fftw3)
