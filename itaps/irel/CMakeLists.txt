set ( MOAB_IREL_SRCS
          ArrayManager.hpp
          AssocPair.cpp
          AssocPair.hpp
          AssocPairSide.hpp
          iRel_Lasso.hpp
          iRel_Lasso.cpp
          Lasso.hpp
          Lasso.cpp
    )
set ( MOAB_IREL_LIB_INCLUDES
          iRel.h
          iRel_f.h
          iRel_protos.h
    )

set (MOAB_DEFINES "")

set(IREL_LINK_LIBS )
if (ENABLE_CGM)
  include_directories( ${CGM_INCLUDES} )
  set (MOAB_DEFINES "-DENABLE_IGEOM ${MOAB_DEFINES}")
  set(IREL_LINK_LIBS ${IGEOM_LIBS} ${CGM_LIBRARIES})

  set(MOAB_IREL_SRCS ${MOAB_IREL_SRCS}
                    GeomAssocPairSide.cpp
                    GeomAssocPairSide.hpp
      )
endif()

set(IREL_LINK_LIBS MOAB ${IREL_LINK_LIBS})
if (ENABLE_FBIGEOM)
  include_directories( ${CMAKE_SOURCE_DIR}/itaps/fbigeom ${CMAKE_BINARY_DIR}/itaps/fbigeom )
  set (MOAB_DEFINES "-DENABLE_FBIGEOM ${MOAB_DEFINES}")
  set(IREL_LINK_LIBS FBiGeomMOAB ${IREL_LINK_LIBS})
  set(MOAB_IREL_SRCS ${MOAB_IREL_SRCS}
                     FBGeomAssocPairSide.cpp
                     FBGeomAssocPairSide.hpp
      )
endif()

if (ENABLE_IMESH)
  include_directories( ${CMAKE_SOURCE_DIR}/itaps/imesh ${CMAKE_BINARY_DIR}/itaps/imesh )
  set (MOAB_DEFINES "-DENABLE_IMESH ${MOAB_DEFINES}")
  set(IREL_LINK_LIBS iMesh ${IREL_LINK_LIBS})
  set(MOAB_IREL_SRCS ${MOAB_IREL_SRCS}
                     MeshAssocPairSide.cpp
                     MeshAssocPairSide.hpp
    )
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/itaps
  ${CMAKE_BINARY_DIR}/itaps
  )

if ( MOAB_HAVE_MPI )
  include_directories( ${CMAKE_SOURCE_DIR}/src/parallel
                       ${CMAKE_BINARY_DIR}/src/parallel)
endif ()

if ( MOAB_HAVE_HDF5 AND HDF5_FOUND )
  include_directories( ${HDF5_INCLUDE_DIRECTORIES}
                       ${CMAKE_SOURCE_DIR}/src/io/mhdf/include )
endif( )

set_source_files_properties(
  ${MOAB_IREL_SRCS} COMPILE_FLAGS "${MOAB_DEFINES}"
  )

add_library(
  iRel
  ${MOAB_IREL_SRCS}
  ${MOAB_IREL_LIB_INCLUDES}
  )

target_link_libraries( iRel ${IREL_LINK_LIBS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )

set_target_properties(
  iRel
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )

install(
  TARGETS iRel
  EXPORT MOABTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
set( MOAB_INSTALL_TARGETS ${MOAB_INSTALL_TARGETS} iRel)

install(
  FILES ${MOAB_IREL_LIB_INCLUDES}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

