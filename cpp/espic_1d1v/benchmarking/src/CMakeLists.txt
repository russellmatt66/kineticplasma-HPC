cmake_minimum_required (VERSION 3.22)
project (espic_1d1v C CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(espic1d1v-bench bench-main.cpp)
target_link_libraries(espic1d1v-bench PRIVATE Eigen3::Eigen)
target_compile_options(espic1d1v-bench PRIVATE -O3 -Wall)