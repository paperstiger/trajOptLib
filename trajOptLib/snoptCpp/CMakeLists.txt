cmake_minimum_required(VERSION 2.8)
project(snoptWrapper)

option(DEBUG "Debug" ON)
if(DEBUG)
    add_definitions(-DDEBUG)
endif(DEBUG)

include_directories(${SNPINC})
include_directories(${PROJECT_SOURCE_DIR}/include)
link_directories(${LIBDIR})
link_directories(${SNPLIBPATH})
link_directories(/usr/lib)
link_directories(/usr/local/lib)
SET(CMAKE_INSTALL_RPATH ${SNPLIBPATH})

include_directories(${TIGER_DIR})
include_directories("/usr/local/include/python2.7")
include_directories("/usr/include/python2.7")
include_directories("/usr/local/include/eigen3")
include_directories("/usr/include/eigen3")
set(server_FILES
    src/toyfunction.cpp src/solverwrapper.cpp src/funcStyle.cpp src/classStyle.cpp
    include/snoptWrapper.h include/funcStyle.h include/classStyle.h)
add_library(snopt SHARED ${server_FILES})
target_link_libraries(snopt python2.7 ${SNPLIB})
set_target_properties(snopt PROPERTIES
       LIBRARY_OUTPUT_DIRECTORY  "${PYLIBDIR}")
# set_property(TARGET snopt PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)


# install(TARGETS snopt DESTINATION ${PYLIBDIR})
