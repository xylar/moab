SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake )

# make sure we have the essential elements installed
FIND_PROGRAM(PYTHON "python")

if(NOT PYTHON)
  FIND_PROGRAM(PYTHON "python3")
endif()

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

include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/pymoab
    ${CMAKE_CURRENT_BINARY_DIR}/pymoab
)

# if all is well go ahead and build/install pymoab interface
if (PYTHON AND NOT PYTHON_MODULES_MISSING)

  ADD_SUBDIRECTORY(pymoab)
  ADD_SUBDIRECTORY(tests)

  #assume a standard install location to start
  SET(PYMOAB_RELATIVE_INSTALL "/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/")
  SET(PYMOAB_USER_INSTALL_LOCATION YES)

if(NOT PYMOAB_PREFIX)
  if(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(PYMOAB_PREFIX "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    #set custom install indicator to true
    SET(PYMOAB_USER_INSTALL_LOCATION YES)
  else()
    # if no prefix is supplied, install to system location
    SET(PYMOAB_PREFIX "/usr/lib")
    SET(PYMOAB_RELATIVE_INSTALL "/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/dist-packages/")
    SET(PYMOAB_USER_INSTALL_LOCATION NO)
  endif(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
endif(NOT PYMOAB_PREFIX)

set(PYMOAB_INSTALL_PREFIX "${PYMOAB_PREFIX}/${PYMOAB_RELATIVE_INSTALL}" CACHE FILEPATH "Location of PyMOAB installation" FORCE)

if(PYMOAB_USER_INSTALL_LOCATION)
      MESSAGE(WARNING "\
The PyMOAB interface has been installed to a custom location. \
In order to access this module from a python intpreter, please add \
its install location to the PYTHONPATH environment variable: \
PYTHONPATH=${PYMOAB_INSTALL_PREFIX}/:\$PYTHONPATH \
")
endif()

  # configure the setup script
  SET(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.cmake.in")
  SET(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY})

  # setup dependencies of the build command
  SET(DEPS   "${CMAKE_CURRENT_SOURCE_DIR}/pymoab/__init__.py")
  SET(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build")

  # python setuptools build/install local version of pymoab for testing
  ADD_CUSTOM_TARGET(pymoab-local-install ALL
    COMMAND export CFLAGS='-w' &&
    ${PYTHON} setup.py build_ext --inplace
    DEPENDS ${DEPS} MOAB)

  # move pymoab egg into the correct location for a --user setuptools install based on username
  SET(INSTALL_SCRIPT_IN "${CMAKE_CURRENT_SOURCE_DIR}/install.sh.cmake")
  SET(INSTALL_SCRIPT    "${CMAKE_CURRENT_BINARY_DIR}/install.sh")
  CONFIGURE_FILE(${INSTALL_SCRIPT_IN} ${INSTALL_SCRIPT})
  INSTALL(CODE "EXECUTE_PROCESS(COMMAND bash ${INSTALL_SCRIPT})")
#  INSTALL(DIRECTORY ${PYMOAB_PREINSTALL_INSTALL_PATH} DESTINATION ${PYMOAB_INSTALL_PREFIX} )

endif()
