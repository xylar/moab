

# create symlinks for all fiels in directory
FILE(GLOB PYMOAB_FILES ${CMAKE_CURRENT_SOURCE_DIR} "*")
IF(NOT ${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_BINARY_DIR})
  FOREACH(PYMOAB_FILE ${PYMOAB_FILES})
    GET_FILENAME_COMPONENT(FILE ${PYMOAB_FILE} NAME)
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${PYMOAB_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${FILE})
  ENDFOREACH(PYMOAB_FILE)
ENDIF()

# add pymoab build target, depends on MOAB target so that all headers are in place before attempting the build
if(ENABLE_DEVMODE)
  include (UseCython)
  cython_add_module(core ${CMAKE_CURRENT_SOURCE_DIR}/core.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
  cython_add_module(hcoord ${CMAKE_CURRENT_SOURCE_DIR}/hcoord.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
  cython_add_module(rng ${CMAKE_CURRENT_SOURCE_DIR}/rng.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
  cython_add_module(scd ${CMAKE_CURRENT_SOURCE_DIR}/scd.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
  cython_add_module(skinner ${CMAKE_CURRENT_SOURCE_DIR}/skinner.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
  cython_add_module(tag ${CMAKE_CURRENT_SOURCE_DIR}/tag.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
  cython_add_module(topo_util ${CMAKE_CURRENT_SOURCE_DIR}/topo_util.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
  cython_add_module(types ${CMAKE_CURRENT_SOURCE_DIR}/types.pyx ${CMAKE_CURRENT_SOURCE_DIR}/core.pxd)
endif(ENABLE_DEVMODE)

