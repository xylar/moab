set(IREL_SOURCES
	AssocPair.cpp
	iRel_Lasso.cpp
	Lasso.cpp
)
set(IREL_INCLUDE_HEADERS
	ArrayManager.hpp
	AssocPair.hpp
	AssocPairSide.hpp
	iRel_Lasso.hpp
	Lasso.hpp
)

add_library(
  iRelMOAB
  ${IREL_SOURCES}
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
  ${CGM_INCLUDE_DIRS}
 )

if ( MOAB_HAVE_MPI )
  include_directories( ${CMAKE_SOURCE_DIR}/src/parallel
                       ${CMAKE_BINARY_DIR}/src/parallel)
endif ()

target_link_libraries( iRelMOAB iMesh MOAB ${IGEOM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )

set_target_properties(
  iRelMOAB
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )
