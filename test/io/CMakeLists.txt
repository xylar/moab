include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/test
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/parallel
)

set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest")

set( TESTS readutil_test.cpp VtkTest.cpp cub_file_test.cpp gmsh_test.cpp nastran_test.cpp ideas_test.cpp stl_test.cpp smf_test.cpp tqdcfr_test.cpp readrtt_test.cpp )
if( MOAB_HAVE_NETCDF )
  set( TESTS ${TESTS} exodus_test.cpp )
  set( TESTS ${TESTS} read_nc.cpp )
  set( TESTS ${TESTS} read_ucd_nc.cpp )
  set( TESTS ${TESTS} read_mpas_nc.cpp )
  set( TESTS ${TESTS} write_nc.cpp )
endif()

if(MOAB_HAVE_CCMIO)
  set( TESTS ${TESTS} ccmio_test.cpp )
endif()

if( MOAB_HAVE_CGNS )
  set( TESTS ${TESTS} cgns_test.cpp )
endif()

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" base ${fname} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES}" )
  target_link_libraries( ${base} MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

