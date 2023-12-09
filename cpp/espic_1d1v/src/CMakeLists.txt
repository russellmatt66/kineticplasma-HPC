cmake_minimum_required(VERSION 3.22)
project(espic1d1v C CXX)

set(CMAKE_CXX_STANDARD 17)

set(EIGEN3_INCLUDE_DIRS "/home/matt/software/eigen-3.4.0")

include_directories(${EIGEN3_INCLUDE_DIRS})

add_executable(espic1d1v main.cpp)
target_compile_options(espic1d1v PRIVATE -O3 -Wall)


