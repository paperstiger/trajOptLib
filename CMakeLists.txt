cmake_minimum_required(VERSION 3.0)
project(trajOptLib)

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
set(TIGER_DIR ${CMAKE_SOURCE_DIR}/lib/TigerTools)

# add_subdirectory(lib/solver_interface)
add_subdirectory(test/rotorWrapper)
