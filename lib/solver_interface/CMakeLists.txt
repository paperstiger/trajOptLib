cmake_minimum_required(VERSION 2.8.3)
project(snopt_cpp)

set(CMAKE_CXX_STANDARD 11)
find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)

set(SNOPT_PATH "/usr/local" CACHE STRING "Where is snopt installed?")
set(IPOPT_PATH "/usr/local" CACHE STRING "Where is snopt installed?")
message("Try to find SNOPT in ${SNOPT_PATH}")
message("Try to find IPOPT in ${IPOPT_PATH}")

set(COMMON_HEADER ${CMAKE_SOURCE_DIR}/common)
include_directories(${EIGEN3_INCLUDE_DIR} ${SNOPT_PATH}/include ${IPOPT_PATH}/include ${COMMON_HEADER})
include_directories(snopt-interface/include ipopt-interface/include)

set(TARGET_LIB_DIRECTORY ${CMAKE_SOURCE_DIR})

link_directories(${SNOPT_PATH}/lib)
link_directories(${IPOPT_PATH}/lib)


set(SNOPTWRAPPERSRC snopt-interface/src/funcStyle.cpp snopt-interface/src/toyfunction.cpp)

set(IPOPTWRAPPERSRC ipopt-interface/src/ipoptWrapper.cpp)

pybind11_add_module(pyoptsolver MODULE ${SNOPTWRAPPERSRC} ${IPOPTWRAPPERSRC} wrapper.cpp)
target_compile_definitions(pyoptsolver PRIVATE SNOPT IPOPT)
if(CMAKE_BUILD_TYPE EQUAL Debug)
    target_compile_definitions(pyoptsolver PRIVATE -DDEBUG)
endif()

target_link_libraries(pyoptsolver PRIVATE snopt7_cpp snopt7 ipopt)
set_target_properties(pyoptsolver PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY  "${TARGET_LIB_DIRECTORY}"
    )
