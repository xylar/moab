cmake_minimum_required( VERSION 2.8.12 )
enable_language(Fortran)

set ( MOAB_DIR "" CACHE PATH "Path to MOAB library installation" )

find_package(MOAB REQUIRED
  HINTS ${MOAB_DIR} ../../ $ENV{MOAB_DIR})

set ( CMAKE_C_COMPILER ${MOAB_CC} )
set ( CMAKE_CXX_COMPILER ${MOAB_CXX} )
set ( CMAKE_Fortran_COMPILER ${MOAB_FC} )

project( MOAB_Examples )

set ( BEGINNER_EXAMPLES
  beginner/mbex1.cpp
  beginner/mbex2.cpp
  beginner/mbex3.cpp
  beginner/mbex4.cpp
)

set ( BASIC_EXAMPLES
  basic/ErrorHandlingSimulation.cpp
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

set ( INTERMEDIATE_EXAMPLES
  intermediate/DirectAccessNoHoles.cpp
  intermediate/DirectAccessWithHoles.cpp
  intermediate/LoadPartial.cpp
  intermediate/PointInElementSearch.cpp
  intermediate/ReduceExchangeTags.cpp
)

set ( ADVANCED_EXAMPLES
  advanced/CrystalRouterExample.cpp
  advanced/ExtrudePoly.cpp
  advanced/GenLargeMesh.cpp
  advanced/LloydRelaxation.cpp
  advanced/VisTags.cpp
)

set ( FORTRAN_EXAMPLES
  fortran/DirectAccessNoHolesF90.F90
  fortran/PushParMeshIntoMoabF90.F90
)

set ( APPS_EXAMPLES
  apps/earth/ContinentsOnGlobe.cpp
)

include_directories(${MOAB_INCLUDE_DIRS} ${MOAB_PACKAGE_INCLUDES})

set ( ALLCXXEXAMPLES ${BEGINNER_EXAMPLES} ${BASIC_EXAMPLES} ${INTERMEDIATE_EXAMPLES} ${ADVANCED_EXAMPLES} ${APPS_EXAMPLES} )

set ( CMAKE_LINKER ${MOAB_CXX} )
foreach( fname ${ALLCXXEXAMPLES} )
  string( REPLACE ".cpp" "" base ${fname} )
  get_filename_component(base ${base} NAME_WE)
  message("-- Driver [${base}] (C++ sources) :: ${fname}")
  add_executable( ${base} ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "'-DMESH_DIR=\"${MOAB_MESH_DIR}\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES}" )
  target_link_libraries( ${base} LINK_PUBLIC ${MOAB_LIBRARIES} )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

#set ( CMAKE_Fortran_COMPILER ${MOAB_FC} )
if ( ENABLE_IMESH )
  foreach( fname ${FORTRAN_EXAMPLES} )
    string( REPLACE ".F90" "" base ${fname} )
    get_filename_component(base ${base} NAME_WE)
    message("-- Driver [${base}] (F90 sources) :: ${fname}")
    add_executable( ${base} ${fname})
    set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "'-DMESH_DIR=\"${MOAB_MESH_DIR}\"' ${TEST_COMP_FLAGS} ${MOAB_DEFINES}" LINKER_LANGUAGE Fortran )
    target_link_libraries( ${base} LINK_PUBLIC ${IMESH_LIBRARIES} )
    add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  endforeach()
endif ( ENABLE_IMESH )

