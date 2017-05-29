# Check that CGM has been enabled

if ( ENABLE_CGM AND MOAB_HAVE_CGM MATCHES ON)
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src )

if ( CGM_FOUND )
  include_directories( ${CGM_INCLUDE_DIRS} )
endif( )


if ( MOAB_HAVE_CGM )
  if ( MOAB_ENABLE_TESTING )
    target_link_libraries( test_geom ${CGM_LIBRARIES} )
  endif()
  option ( MOAB_BUILD_CGM2MOAB    "Build cgm2moab?"  OFF )
  if(MOAB_BUILD_CGM2MOAB)
    add_executable( cgm2moab cgm2moab.cpp
                             obb_analysis.cpp
                             cgm2moab.hpp )
    add_executable( cub2h5m  cub2h5m.cpp
    		             quads_to_tris.cpp
			     quads_to_tris.hpp)
    target_link_libraries( cgm2moab
                           MOAB
                           ${CGM_LIBRARIES} )
    target_link_libraries( cub2h5m
                           MOAB
                           ${CGM_LIBRARIES} )
    install(
      TARGETS cgm2moab
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )
    install(
      TARGETS cub2h5m
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )
  endif()
endif( )

endif ( ENABLE_CGM AND MOAB_HAVE_CGM MATCHES ON)
