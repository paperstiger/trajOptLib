cmake_minimum_required(VERSION 3.0)
project(snoptWrapper)

option(DEBUG "Debug" OFF)
if(DEBUG)
    add_definitions(-DDEBUG)
endif(DEBUG)

include_directories(${SNPINC})
include_directories(${PROJECT_SOURCE_DIR}/include)
link_directories(${LIBDIR})
link_directories(${SNPLIBPATH})
SET(CMAKE_INSTALL_RPATH ${SNPLIBPATH})

include_directories(${TIGER_DIR})
message("Python version is ${PYBIND11_PYTHON_VERSION}")

set(server_FILES
    src/toyfunction.cpp src/solverwrapper.cpp src/funcStyle.cpp src/classStyle.cpp
    include/snoptWrapper.h include/funcStyle.h include/classStyle.h)
pybind11_add_module(snopt MODULE ${server_FILES})
target_link_libraries(snopt PRIVATE ${SNPLIB})

include_directories(${EIGEN_DIR})

set_target_properties(snopt PROPERTIES
       LIBRARY_OUTPUT_DIRECTORY  "${PYLIBDIR}"
       PREFIX "lib")
