find_program(PYTHON "python")

if (PYTHON)
  SET(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  SET(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  SET(MOAB_PXD_IN "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/moab.pxd.in")
  SET(MOAB_PXD "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/moab.pxd") 
#  SET(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py")
  SET(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build")

  configure_file(${SETUP_PY_IN} ${SETUP_PY})
  configure_file(${MOAB_PXD_IN} ${MOAB_PXD})
  
  add_custom_command(OUTPUT ${OUTPUT}
                     COMMAND ${PYTHON}
                     ARGS setup.py build
		     DEPENDS ${DEPS})
  add_custom_target(target ALL DEPENDS ${OUTPUT})

  install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --user)")
endif()
