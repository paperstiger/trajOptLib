cmake_minimum_required(VERSION 2.8.3)
project(snopt_cpp)

set(CMAKE_CXX_STANDARD 11)
find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)

option(BUILD_SNOPT "Build with SNOPT" ON)
option(BUILD_IPOPT "Build with IPOPT" ON)

set(SRC wrapper.cpp)
set(LINKLIB )

if(BUILD_SNOPT)
    set(SNOPT_PATH "/usr/local" CACHE STRING "Where is snopt installed?")
    message("Try to find SNOPT in ${SNOPT_PATH}")
    include_directories(${SNOPT_PATH}/include)
    link_directories(${SNOPT_PATH}/lib)
    set(LINKLIB ${LINKLIB} snopt7_cpp snopt7)
    set(SNOPTWRAPPERSRC snopt-interface/src/funcStyle.cpp snopt-interface/src/toyfunction.cpp)
    set(SRC ${SRC} ${SNOPTWRAPPERSRC})
endif(BUILD_SNOPT)
if(BUILD_IPOPT)
    set(IPOPT_PATH "/usr/local" CACHE STRING "Where is snopt installed?")
    message("Try to find IPOPT in ${IPOPT_PATH}")
    include_directories(${IPOPT_PATH}/include)
    link_directories(${IPOPT_PATH}/lib)
    set(LINKLIB ${LINKLIB} libipopt.so)
    set(IPOPTWRAPPERSRC ipopt-interface/src/ipoptWrapper.cpp)
    set(SRC ${SRC} ${IPOPTWRAPPERSRC})
endif(BUILD_IPOPT)

if(NOT (BUILD_SNOPT OR BUILD_IPOPT))
    message(FATAL_ERROR "You have to build at least one of SNOPT and IPOPT")
endif()

set(COMMON_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/common)
include_directories(${EIGEN3_INCLUDE_DIR} ${COMMON_HEADER})
include_directories(snopt-interface/include ipopt-interface/include)

set(TARGET_LIB_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

pybind11_add_module(pyoptsolvercpp MODULE ${SRC})
if(BUILD_SNOPT)
    target_compile_definitions(pyoptsolvercpp PRIVATE SNOPT)
endif()
if(BUILD_IPOPT)
    target_compile_definitions(pyoptsolvercpp PRIVATE IPOPT ENABLEIP)
    target_compile_definitions(pyoptsolvercpp PRIVATE HAVE_CSTDDEF)
endif()

if(CMAKE_BUILD_TYPE EQUAL Debug)
    target_compile_definitions(pyoptsolvercpp PRIVATE -DDEBUG)
endif()

target_link_libraries(pyoptsolvercpp PRIVATE ${LINKLIB})
set_target_properties(pyoptsolvercpp PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY  "${TARGET_LIB_DIRECTORY}"
    )
