####################################################################################
#       MOAB - CMake downstream application configuration/compilation template     #
#                                                                                  #
#  Configuration Usage:                                                            #
#       cmake -DMOAB_ROOT=${PATH_TO_MOAB_INSTALLATION} ${EXAMPLE_SRC_PATH}         #
####################################################################################
##                           DO NOT MODIFY BELOW                                  ##
####################################################################################
cmake_minimum_required ( VERSION 2.8.12 )

set ( MOAB_ROOT "" CACHE PATH "Path to MOAB library installation" )

find_package(MOAB REQUIRED
  HINTS ${MOAB_ROOT} ${MOAB_ROOT}/lib/cmake/MOAB)

#Add our Cmake directory to the module search path
set ( CMAKE_MODULE_PATH ${MOAB_DIR} ${CMAKE_MODULE_PATH} )

set ( CMAKE_C_COMPILER ${MOAB_CC} )
set ( CMAKE_CXX_COMPILER ${MOAB_CXX} )
set ( CMAKE_Fortran_COMPILER ${MOAB_FC} )
set ( CMAKE_C_FLAGS ${MOAB_C_FLAGS} )
set ( CMAKE_CXX_FLAGS ${MOAB_CXX_FLAGS} )
set ( CMAKE_Fortran_FLAGS ${MOAB_FC_FLAGS} )
set ( CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} )

if ( ENABLE_IMESH )
  enable_language(Fortran)
endif ( ENABLE_IMESH )

message ( "++ Found MOAB configuration details at: ${MOAB_DIR}" )

####################################################################################
##                      MODIFY BELOW AS NECESSARY                                 ##
##     (Replace `ExampleMOABProgram` with your source file name)                  ##
####################################################################################

project (MOAB_CMAKE_TEMPLATE_EXAMPLE)

include_directories(${MOAB_INCLUDE_DIRS})

add_executable(ExampleMOABProgram ExampleMOABProgram.cpp)
set_target_properties( ExampleMOABProgram PROPERTIES COMPILE_FLAGS "'-DMESH_DIR=\"${MOAB_MESH_DIR}\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES}" )
target_link_libraries( ExampleMOABProgram LINK_PUBLIC ${MOAB_LIBRARIES} )

