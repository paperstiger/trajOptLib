cmake_minimum_required (VERSION 2.8)

# projectname is the same as the main-executable
project(solver_interface)

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
if(NOT pybind11_FOUND)
    add_subdirectories(${CMAKE_SOURCE_DIR}/../pybind11)
endif()

add_subdirectory(pyoptsolver/src)
