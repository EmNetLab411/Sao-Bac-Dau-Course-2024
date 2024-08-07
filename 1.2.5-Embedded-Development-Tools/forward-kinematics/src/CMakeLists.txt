cmake_minimum_required(VERSION 3.13)
project(forward_kinematics)
set(CMAKE_CXX_FLAGS "-Wall -fopenmp")
add_executable(out main.cpp)
# Link OpenMP
target_link_libraries(out ${OpenMP_CXX_LIBRARIES})
# Link Google benchmark
find_package(benchmark REQUIRED)
target_link_libraries(out benchmark)
# Link linalg_utils
target_link_directories(out PRIVATE ${CMAKE_SOURCE_DIR}/linalg/build/)
target_link_libraries(out linalg)
