# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (MOAB_EXAMPLES)

set ( MOAB_DIR "" CACHE PATH "Path to MOAB library installation" )
find_package(MOAB REQUIRED)

#Add our Cmake directory to the module search path
set(CMAKE_MODULE_PATH ${MOAB_DIR} ${CMAKE_MODULE_PATH})

include_directories(${MOAB_INCLUDE_DIRS})

add_executable(HelloMOAB HelloMOAB.cpp)
target_link_libraries(HelloMOAB LINK_PUBLIC ${MOAB_LIBRARIES})

