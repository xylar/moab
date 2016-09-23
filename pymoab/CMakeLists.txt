find_program(PYTHON "python")
find_program(NOSE "nosetests")

if (PYTHON AND NOSE)
  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pymoab" FILES_MATCHING PATTERN "*")
  SET(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  SET(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  SET(MOAB_PXD_IN "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/moab.pxd.in")
  SET(MOAB_PXD "${CMAKE_CURRENT_BINARY_DIR}/pymoab/moab.pxd") 
  SET(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/__init__.py")
  SET(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build")

  configure_file(${SETUP_PY_IN} ${SETUP_PY})
  configure_file(${MOAB_PXD_IN} ${MOAB_PXD})
  
  add_custom_command(OUTPUT ${OUTPUT}
                     COMMAND ${PYTHON}
                     ARGS setup.py build
		     DEPENDS ${DEPS})
  add_custom_target(pymoab-build ALL DEPENDS ${OUTPUT})

  add_custom_target(pymoab-local-install ALL COMMAND export PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$PYTHONPATH && ${PYTHON} setup.py install --install-platlib=${CMAKE_CURRENT_BINARY_DIR} DEPENDS pymoab-build)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pymoab-0.0.0-py2.7-linux-x86_64.egg DESTINATION /home/$ENV{USER}/.local/lib/python2.7/site-packages/)

  add_custom_target(pymoab-tests COMMAND cd ${CMAKE_CURRENT_BINARY_DIR}/pymoab/test/ && nosetests)
endif()

