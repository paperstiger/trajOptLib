cmake_minimum_required(VERSION 3.0)
project(rotorWrapper)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${TIGER_DIR})
message("Tiger dir is ${TIGER_DIR}")

set(Rotor_FILES
    src/rotorWrapper.cpp
    )
pybind11_add_module(Rotor MODULE ${Rotor_FILES})

set_target_properties(Rotor PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY  "${CMAKE_SOURCE_DIR}/test"
    PREFIX "lib")

install(TARGETS Rotor
    DESTINATION ${CMAKE_SOURCE_DIR}/test)
