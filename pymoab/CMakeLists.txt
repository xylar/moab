
# make sure we have the essential elements installed
FIND_PROGRAM(PYTHON "python")
FIND_PROGRAM(NOSE "nosetests")

# make sure that specific python modules are installed
SET(PYTHON_MODULES_EXIST false)
EXECUTE_PROCESS(COMMAND "${PYTHON}" "-c" "import numpy" RESULT_VARIABLE PYTHON_MODULES_MISSING)
EXECUTE_PROCESS(COMMAND "${PYTHON}" "-c" "import Cython" RESULT_VARIABLE PYTHON_MODULES_MISSING)
EXECUTE_PROCESS(COMMAND "${PYTHON}" "-c" "import setuptools" RESULT_VARIABLE PYTHON_MODULES_MISSING)

if (PYTHON_MODULES_MISSING)
  MESSAGE( SEND_ERROR "Could not detect necessary python modules.")
endif()

# if all is well go ahead and build/install pymoab interface
if (PYTHON AND NOSE AND NOT PYTHON_MODULES_MISSING)
  # copy all source to current dir so .cpp's aren't built in source
  FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pymoab" FILES_MATCHING PATTERN "*")

  # configure the setup script
  SET(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  SET(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY})

  # configure moab pxd
  SET(MOAB_PXD_IN "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/moab.pxd.in")
  SET(MOAB_PXD "${CMAKE_CURRENT_BINARY_DIR}/pymoab/moab.pxd") 
  CONFIGURE_FILE(${MOAB_PXD_IN} ${MOAB_PXD})

  # setup dependencies of the build command
  SET(DEPS   "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/__init__.py")
  SET(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build")
  # build pymoab
  ADD_CUSTOM_COMMAND(OUTPUT ${OUTPUT}
                     COMMAND ${PYTHON}
                     ARGS setup.py build
		     DEPENDS ${DEPS})
  # add pymoab build target, depends on MOAB target so that all headers are in place before attempting the build
  ADD_CUSTOM_TARGET(pymoab-build ALL DEPENDS ${OUTPUT} MOAB)

  # pymoab install- depends on pymoab-build and is run during make
  # !!!!THIS INSTALL IS DONE LOCALLY AND COPIED TO THE CORRECT LOCATOIN DURING INSTALL !!!! (Subject to Change)
  ADD_CUSTOM_TARGET(pymoab-local-install ALL COMMAND export PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$PYTHONPATH && ${PYTHON} setup.py install --install-platlib=${CMAKE_CURRENT_BINARY_DIR} DEPENDS pymoab-build)
  # move pymoab egg into the correct location for a --user setuptools install based on username
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/pymoab-0.0.0-py2.7-linux-x86_64.egg DESTINATION /home/$ENV{USER}/.local/lib/python2.7/site-packages/)

  # run all pymoab tests, local install is added to PYTHONPATH so that make check can be run before install
  ADD_CUSTOM_TARGET(pymoab-tests COMMAND export PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH} && cd ${CMAKE_CURRENT_BINARY_DIR}/pymoab/test/ && nosetests)
endif()

