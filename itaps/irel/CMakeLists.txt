set( IREL_SOURCE Lasso.cpp iRel_Lasso.cpp AssocPair.cpp FBGeomAssocPairSide.cpp GeomAssocPairSide.cpp )
set(
  IREL_INCLUDE_HEADERS
  Lasso.hpp
  iRel_Lasso.hpp
  iRel.h
  iRel_protos.h
  iRel_f.h
  ArrayManager.hpp
  AssocPair.hpp
  AssocPairSide.hpp
  FBGeomAssocPairSide.hpp
  GeomAssocPairSide.hpp
  MeshAssocPairSide.hpp
  )

add_library(
  iRelMOAB
  ${IREL_SOURCE}
  ${IREL_INCLUDE_HEADERS}
  )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/itaps
  ${CMAKE_BINARY_DIR}/itaps
  ${CMAKE_SOURCE_DIR}/itaps/imesh
  ${CMAKE_BINARY_DIR}/itaps/imesh
  ${CMAKE_SOURCE_DIR}/itaps/fbigeom
  ${CMAKE_BINARY_DIR}/itaps/fbigeom
  ${CGM_INCLUDES}
)

if ( MOAB_HAVE_MPI )
  include_directories( ${CMAKE_SOURCE_DIR}/src/parallel
                       ${CMAKE_BINARY_DIR}/src/parallel)
endif ()

target_link_libraries( iRelMOAB FBiGeomMOAB iMesh MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )

set_target_properties(
  iRelMOAB
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )

