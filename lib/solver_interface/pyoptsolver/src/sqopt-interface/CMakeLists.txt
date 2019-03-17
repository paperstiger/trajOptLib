cmake_minimum_required(VERSION 2.8.3)
project(sqopt_cpp)

set(CMAKE_CXX_STANDARD 11)

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

set(SNOPT_PATH "/usr/local" CACHE STRING "Where is snopt installed?")
set(TARGET_LIB_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(SNOPT_LIB_DIRECTORY ${SNOPT_PATH}/lib)

link_directories(${SNOPT_LIB_DIRECTORY})
include_directories(include)

add_library(sqopt_cpp SHARED src/sqoptEZ.cpp)

add_executable(sqopt_example cppexamples/hs118ez.cpp)
target_link_libraries(sqopt_example sqopt_cpp sqopt7 snopt7_cpp snopt7)

pybind11_add_module(sqopt MODULE src/sqoptEZ.cpp src/sqopt_wrapper.cpp)
target_include_directories(sqopt PRIVATE include ${EIGEN3_INCLUDE_DIR})
target_link_libraries(sqopt PRIVATE sqopt_cpp sqopt7 snopt7_cpp snopt7)
set_target_properties(sqopt PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY  "${TARGET_LIB_DIRECTORY}"
       PREFIX "lib")

install(TARGETS sqopt_cpp sqopt
    DESTINATION ${TARGET_LIB_DIRECTORY}
    )
install(FILES include/sqopt.h include/sqoptProblem.hpp include/sqoptEZ.h
    DESTINATION ${TARGET_LIB_DIRECTORY}/include}
    )
