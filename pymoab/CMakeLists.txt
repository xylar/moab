SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake )

# make sure we have the essential elements installed
FIND_PROGRAM(PYTHON "python")

find_package(PythonInterp REQUIRED)
find_package(PythonDev REQUIRED)
# make sure that specific python modules are installed
SET(NUMPY_MODULE_MISSING false)
SET(SETUPTOOLS_MODULE_MISSING false)

EXECUTE_PROCESS(COMMAND "${PYTHON}" "-c" "import numpy" RESULT_VARIABLE NUMPY_MODULE_MISSING)
if (NUMPY_MODULE_MISSING)
  MESSAGE( SEND_ERROR "Python numpy module is missing.")
endif()
EXECUTE_PROCESS(COMMAND "${PYTHON}" "-c" "import setuptools" RESULT_VARIABLE SETUPTOOLS_MODULE_MISSING)
if (SETUPTOOLS_MODULE_MISSING)
  MESSAGE( SEND_ERROR "Python setuptools module is missing.")
endif()

if(MOAB_FORCE_32_BIT_HANDLES)
  MESSAGE( SEND_ERROR "Cannot build pymoab with 32bit EntityHandles enabled. Please remove the --enable-32bit-handles configuration option or build without MOAB's python interface.")
endif()

if(NOT BUILD_SHARED_LIBS)
  MESSAGE( SEND_ERROR "Cannot build static libraries and the python interface. Please build with shared libraries (-DBUILD_SHARED_LIBS=ON) or build without MOAB's python interface.")
endif()

# if all is well go ahead and build/install pymoab interface
if (PYTHON AND NOT PYTHON_MODULES_MISSING)
  
  ADD_SUBDIRECTORY(pymoab)
  ADD_SUBDIRECTORY(tests)

  #assume a standard install location to start
  SET(PYMOAB_CUSTOM_INSTALL_LOCATION NO)

  if(PYMOAB_PREFIX)
    #if a pymoab prefix has been set, us it
    SET(PYMOAB_INSTALL_PREFIX ${PYMOAB_PREFIX})
    #set custom install indicator to true
    SET(PYMOAB_CUSTOM_INSTALL_LOCATION YES)
  elseif(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    #otherwise install to the same location as the CMAKE_INSTALL_PREFIX
    SET(PYMOAB_INSTALL_PREFIX "${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/")
    #set custom install indicator to true    
    SET(PYMOAB_CUSTOM_INSTALL_LOCATION YES)
  else()
    # if no prefix is supplied, install to system location
    SET(PYMOAB_INSTALL_PREFIX "/usr/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/dist-packages/")
  endif()

  if(PYMOAB_CUSTOM_INSTALL_LOCATION)
      MESSAGE(WARNING "\
The PyMOAB interface has been installed to a custom location. \
In order to access this module from a python intpreter, please add \
its install location to the PYTHONPATH environment variable:
\
PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${PYMOAB_INSTALL_PREFIX}:\$PYTHONPATH \
")
  endif()
    
  # configure the setup script
  SET(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.cmake.in")
  SET(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY})

  # setup dependencies of the build command
  SET(DEPS   "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/__init__.py")
  SET(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build")
  
  # add pymoab build target, depends on MOAB target so that all headers are in place before attempting the build
  # ADD_CUSTOM_TARGET(pymoab-build ALL DEPENDS ${OUTPUT} pymoab-cythonize MOAB)

  # build pymoab
  ADD_CUSTOM_TARGET(pymoab-build ALL COMMAND export CFLAGS='-w' && ${PYTHON} setup.py build DEPENDS ${DEPS} MOAB)

  #create local installation directory
  ADD_CUSTOM_TARGET(pymoab-local-install-dir ALL COMMAND ${CMAKE_COMMAND} -E make_directory ${PYMOAB_LOCAL_INSTALL_PATH} DEPENDS pymoab-build)
  # pymoab install- depends on pymoab-build and is run during make
  # !!!!THIS INSTALL IS DONE LOCALLY AND COPIED TO THE CORRECT LOCATOIN DURING INSTALL !!!! (Subject to Change)
  ADD_CUSTOM_TARGET(pymoab-local-install ALL COMMAND export PYTHONPATH=${PYMOAB_LOCAL_INSTALL_PATH}:$PYTHONPATH && ${PYTHON} setup.py install --install-platlib=${PYMOAB_LOCAL_INSTALL_PATH} --install-purelib=${PYMOAB_LOCAL_INSTALL_PATH} DEPENDS pymoab-local-install-dir)
  # move pymoab egg into the correct location for a --user setuptools install based on username
  INSTALL(DIRECTORY ${PYMOAB_LOCAL_INSTALL_PATH}  DESTINATION ${PYMOAB_INSTALL_PREFIX} )
  
endif()

