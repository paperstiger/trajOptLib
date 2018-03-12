cmake_minimum_required(VERSION 2.8)
project(snoptWrapper)

option(DEBUG "Debug" OFF)
if(DEBUG)
    add_definitions(-DDEBUG)
endif(DEBUG)

include_directories(${SNPINC})
include_directories(${PROJECT_SOURCE_DIR}/include)
link_directories(${LIBDIR})
link_directories(${SNPLIBPATH})
link_directories(/usr/lib)
link_directories(/usr/local/lib)


include_directories("/usr/local/include/python2.7")
include_directories("/usr/include/python2.7")
include_directories("/usr/local/include/eigen3")
include_directories("/usr/include/eigen3")
set(server_FILES
    src/toyfunction.cpp src/solverwrapper.cpp src/funcStyle.cpp src/classStyle.cpp
    include/snoptWrapper.h include/funcStyle.h include/classStyle.h)
add_library(snoptWrapper SHARED ${server_FILES})
target_link_libraries(snoptWrapper python2.7 ${PYLIBDIR})

install(TARGETS snopt
    DESTINATION ${LIBDIR})
