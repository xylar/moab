add_subdirectory(point_location)
set( LIBS MOAB )
set( TESTS seqperf.cpp
           adj_time.cpp
           adj_mem_time_test.cpp
           )

set( TOOLS perf.cpp
           perftool.cpp
           umr_perf.cpp )

if ( MOAB_HAVE_IMESH )
  set(TESTS ${TESTS} tstt_perf_binding.cpp)
  set(LIBS ${LIBS} iMesh )
endif()
if ( MOAB_HAVE_CGM )
  set(LIBS ${LIBS} ${CGM_LIBRARIES})
endif()

set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest")

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" tmp ${fname} )
  string( REPLACE ".cc" "" base ${tmp} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES} ${TEST_COMP_FLAGS} -DTEST" )
  target_link_libraries( ${base} ${LIBS} )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

foreach( fname ${TOOLS} )
  string( REPLACE ".cpp" "" tmp ${fname} )
  string( REPLACE ".cc" "" base ${tmp} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES} ${TEST_COMP_FLAGS} -DTEST" )
  target_link_libraries( ${base} ${LIBS} )
endforeach()
