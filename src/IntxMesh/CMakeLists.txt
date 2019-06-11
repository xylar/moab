if( MOAB_HAVE_HDF5 )
  include_directories(${HDF5_INCLUDES})
  set( MOAB_HDF_LIBS ${HDF5_LIBRARIES} )
endif()
if( MOAB_HAVE_NETCDF )
  include_directories(${NETCDF_INCLUDES})
  if ( MOAB_HAVE_MPI )
    set( MOAB_NETCDF_LIBS ${PNETCDF_LIBRARIES} ${NETCDF_LIBRARIES} )
  else ( MOAB_HAVE_MPI )
    set( MOAB_NETCDF_LIBS ${NETCDF_LIBRARIES} )
  endif ( MOAB_HAVE_MPI )
endif()

set(MOAB_INTX_LIB_SRCS Intx2Mesh.cpp Intx2MeshInPlane.cpp Intx2MeshOnSphere.cpp IntxRllCssphere.cpp IntxUtils.cpp )
#if(ENABLE_IMESH AND MOAB_HAVE_MPI AND MPI_FOUND)
# include_directories( ${CMAKE_SOURCE_DIR}/itaps/imesh  ${CMAKE_SOURCE_DIR}/itaps )
# set(MOAB_INTX_LIB_SRCS ${MOAB_INTX_LIB_SRCS} wrap_intx.cpp )
#endif()

add_library( moab_intx_mesh OBJECT ${MOAB_INTX_LIB_SRCS} ${MOAB_INTX_INC_HDRS} )
set_target_properties( moab_intx_mesh
                       PROPERTIES COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}")

include(AutoconfHeader)
moab_install_headers(${MOAB_INTX_INC_HDRS})

