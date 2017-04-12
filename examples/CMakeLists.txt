cmake_minimum_required( VERSION 2.8.12 )

project( MOAB_Examples )

set ( MOAB_DIR "" CACHE PATH "Path to MOAB library installation" )

find_package(MOAB REQUIRED)

set ( BEGINNER_EXAMPLES beginner/mbex1.cpp
  beginner/mbex2.cpp
  beginner/mbex3.cpp
  beginner/mbex4.cpp
)

set ( BASIC_EXAMPLES basic/ErrorHandlingSimulation.cpp
  basic/GetEntities.cpp
  basic/HelloMOAB.cpp
  basic/HelloParMOAB.cpp
  basic/ReadWriteTest.cpp
  basic/SetsNTags.cpp
  basic/StructuredMeshSimple.cpp
  basic/TestErrorHandling.cpp
  basic/TestErrorHandlingPar.cpp
  basic/TestExodusII.cpp
  basic/UniformRefinement.cpp
)

include_directories(${MOAB_INCLUDE_DIRS} ${MOAB_PACKAGE_INCLUDES})
#add_executable(HelloMOAB HelloMOAB.cpp)
#target_link_libraries(HelloMOAB LINK_PUBLIC ${MOAB_LIBRARIES})

set ( ALLEXAMPLES ${BEGINNER_EXAMPLES} ${BASIC_EXAMPLES} )

foreach( fname ${ALLEXAMPLES} )
  string( REPLACE ".cpp" "" base ${fname} )
  get_filename_component(base ${base} NAME_WE)
  message("-- Driver [${base}] :: source :: ${fname}")
  add_executable( ${base} ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "'-DMESH_DIR=\"/Users/vijaysm/code/sigma/moab-pristine/MeshFiles/unittest\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES}" )
  target_link_libraries( ${base} LINK_PUBLIC ${MOAB_LIBRARIES} )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()


