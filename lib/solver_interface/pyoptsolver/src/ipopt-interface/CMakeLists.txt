cmake_minimum_required(VERSION 2.8.3)
project(ipopt_cpp)

set(CMAKE_CXX_STANDARD 11)
find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)

set(PYBIND11_PYTHON_VERSION 2.7 CACHE STRING "What is python version")
set(COMMON_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/../common)

include_directories(include ${EIGEN3_INCLUDE_DIR} ${COMMON_HEADER})

set(TARGET_LIB_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)

file(GLOB IPOPTWRAPPERSRC src/*.cpp)

pybind11_add_module(pyipopt MODULE ${IPOPTWRAPPERSRC})
if(CMAKE_BUILD_TYPE EQUAL Debug)
    target_compile_definitions(pyipopt PRIVATE -DDEBUG)
endif()
target_compile_definitions(pyipopt PRIVATE HAVE_CSTDDEF)

target_link_libraries(pyipopt PRIVATE ipopt)
set_target_properties(pyipopt PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY  "${TARGET_LIB_DIRECTORY}"
    PREFIX "lib")
