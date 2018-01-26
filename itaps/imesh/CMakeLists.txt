
set ( MOAB_IMESH_SRCS
          iMesh_MOAB.hpp iMesh_MOAB.cpp
          MBIter.hpp )
set ( MOAB_IMESH_LIB_INCLUDES
          iMesh.h
          iMesh_extensions.h
          iMesh_f.h
          iMesh_protos.h
          MBiMesh.hpp )

if ( MOAB_HAVE_MPI )
  include_directories( ${CMAKE_SOURCE_DIR}/src/parallel
                       ${CMAKE_BINARY_DIR}/src/parallel)
  set(MOAB_IMESH_SRCS ${MOAB_IMESH_SRCS}
                      iMeshP_MOAB.cpp)
  set(MOAB_IMESH_LIB_INCLUDES
          ${MOAB_IMESH_LIB_INCLUDES}
          iMeshP.h
          iMeshP_f.h
          iMeshP_protos.h )
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/itaps
  ${CMAKE_BINARY_DIR}/itaps
  )

if ( MOAB_HAVE_HDF5 AND HDF5_FOUND )
  include_directories( ${HDF5_INCLUDE_DIRECTORIES}
                       ${CMAKE_SOURCE_DIR}/src/io/mhdf/include )
endif( )

set_source_files_properties(
  ${MOAB_IMESH_SRCS} COMPILE_FLAGS "${MOAB_DEFINES}"
  )

add_library(
  iMesh
  ${MOAB_IMESH_SRCS}
  ${MOAB_IMESH_LIB_INCLUDES}
  )

target_include_directories(iMesh
  PUBLIC
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

target_link_libraries( iMesh MOAB ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )

set_target_properties(
  iMesh
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )

install(
  TARGETS iMesh
  EXPORT MOABTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
set( MOAB_INSTALL_TARGETS ${MOAB_INSTALL_TARGETS} iMesh)

install(
  FILES ${MOAB_IMESH_LIB_INCLUDES}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
