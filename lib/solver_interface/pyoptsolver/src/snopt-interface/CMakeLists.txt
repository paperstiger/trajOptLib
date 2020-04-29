cmake_minimum_required(VERSION 2.8.3)
project(snopt_cpp)

set(CMAKE_CXX_STANDARD 11)
find_package(Eigen3 REQUIRED)
set(SNOPT_PATH "/usr/local" CACHE STRING "Where is snopt installed?")
message("Try to find SNOPT in ${SNOPT_PATH}")
include_directories(include ${EIGEN3_INCLUDE_DIR} ${SNOPT_PATH}/include)

set(TARGET_LIB_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
link_directories(${SNOPT_PATH}/lib)


file(GLOB SNOPTWRAPPERSRC src/*.cpp)

# add_executable(sqopt_example cppexamples/hs118ez.cpp)
# target_link_libraries(sqopt_example sqopt_cpp sqopt7 snopt7_cpp snopt7)

# build binding for SNOPT
find_package(pybind11 REQUIRED)

pybind11_add_module(pysnopt MODULE ${SNOPTWRAPPERSRC})
if(CMAKE_BUILD_TYPE EQUAL Debug)
    target_compile_definitions(pysnopt PRIVATE -DDEBUG)
endif()

target_link_libraries(pysnopt PRIVATE snopt7_cpp snopt7)
set_target_properties(pysnopt PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY  "${TARGET_LIB_DIRECTORY}"
    PREFIX "lib")

#install(FILES include/snoptProblem.hh
#    DESTINATION ${TARGET_LIB_DIRECTORY}/include}
#    )
