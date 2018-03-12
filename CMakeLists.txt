cmake_minimum_required(VERSION 2.8)
project(trajOptLib)

# set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_STANDARD 11)

option(DEBUG "Debug" ON)
if(DEBUG)
    add_definitions(-DDEBUG)
    option(CMAKE_BUILD_TYPE "BType" Debug)
endif(DEBUG)

# prepare two libs, snopt and f2c

set(TIGER_DIR "${PROJECT_SOURCE_DIR}/lib")  # so I can #include “TigerTools/TigerEigen.h”
set(SNOPTDIR "snopt7")
set(SNOPTPATH "${PROJECT_SOURCE_DIR}/lib/${SNOPTDIR}")
set(SNPLIB "snoptc")
set(PYLIBDIR "${PROJECT_SOURCE_DIR}/trajOptLib")
set(f2cPATH "${SNOPTPATH}/win64")  # this might change
set(f2cLIB f2c)
# set(INCDIR "${PROJECT_SOURCE_DIR}/Headers")
set(LIBDIR "${PROJECT_SOURCE_DIR}/lib")
# Build a SNOPT library if necessary
file(GLOB snoptcsrc
	"${SNOPTPATH}/csrc/*.c"
	"${SNOPTPATH}/cppsrc/*.h"
	"${SNOPTPATH}/cppsrc/*.hh"
	"${SNOPTPATH}/cppsrc/*.c"
	"${SNOPTPATH}/cppsrc/*.cpp"
	"${SNOPTPATH}/cppsrc/*.cc"
	)
add_library(${SNPLIB}
		SHARED
		${snoptcsrc}
	)
set_target_properties(${SNPLIB} PROPERTIES
			 LIBRARY_OUTPUT_DIRECTORY  "${LIBDIR}")
target_include_directories(${SNPLIB} PUBLIC "${f2cPATH}")
target_link_libraries(${SNPLIB} f2c)
set(SNPLIBPATH ${LIBDIR})
set(SNPINC ${SNOPTPATH}/cppsrc)

add_subdirectory(lib/libf2c)
add_subdirectory(trajOptLib/snoptCpp)  # so I can make libsnopt which contains class definition


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_BUILD_DIRECTORY ${PROJECT_SOURCE_DIR}/Build)

include_directories(${SNPINC})
link_directories(${LIBDIR})
link_directories(${SNPLIBPATH})  # although it is the same with previous one


