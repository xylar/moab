cmake_minimum_required( VERSION 3.0 )
SET(PACKAGE_NAME "MOAB")
SET(PACKAGE_VERSION "5.1.1")

if (NOT WIN32)
  #This all breaks on windows.
  SET(CMAKE_Fortran_COMPILER_INIT ${CMAKE_GENERATOR_FC})
  SET(CMAKE_CXX_FLAGS_INIT "-fPIC -DPIC")
  SET(CMAKE_CXX_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  SET(CMAKE_C_FLAGS_INIT "-fPIC -DPIC")
  SET(CMAKE_C_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_C_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  SET(CMAKE_Fortran_FLAGS_INIT "-fPIC")
  SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O2")
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS_INIT "")
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS_INIT "")
endif()

project( MOAB )

include(GNUInstallDirs)

#Add our Cmake directory to the module search path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/config ${CMAKE_MODULE_PATH})
################################################################################
# Set up version info
################################################################################
include (config/GetAcInitVersion.cmake)
get_ac_init_version()
set ( MOAB_VERSION_STRING "${PACKAGE_NAME} ${VERSION_STRING}" )
set ( MOAB_VERSION        "${VERSION_STRING}" )
set ( MOAB_VERSION_MAJOR  ${MAJOR_VERSION}  )
set ( MOAB_VERSION_MINOR  ${MINOR_VERSION}  )
if ( DEFINED PATCH_VERSION )
  set ( MOAB_VERSION_PATCH "${PATCH_VERSION}" )
else ( DEFINED PATCH_VERSION )
  if ( MOAB_VERSION_MINOR EQUAL 99 )
    set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (alpha)" )
  else ( MOAB_VERSION_MINOR EQUAL 99 )
    set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (beta)" )
  endif ( MOAB_VERSION_MINOR EQUAL 99 )
endif ( DEFINED PATCH_VERSION )

IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR CMAKE_VERSION VERSION_GREATER "3.0.0")
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0020 NEW)
  cmake_policy(SET CMP0042 OLD)
ENDIF ()

set ( abs_srcdir   ${CMAKE_SOURCE_DIR} )
set ( abs_builddir ${PROJECT_BINARY_DIR} )

################################################################################
# Install Related Settings
################################################################################

get_filename_component( MOAB_ABSSRC_DIR2 moab.make.in REALPATH )
execute_process( COMMAND dirname ${MOAB_ABSSRC_DIR2}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE MOAB_ABSSRC_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

## Set the directory where the binaries will be stored
set( EXECUTABLE_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Directory where all executable will be stored"
)

## Set the directory where the libraries will be stored
set( LIBRARY_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH
  "Directory where all the libraries will be stored"
)
mark_as_advanced(
  MOAB_ABSSRC_DIR
  EXECUTABLE_OUTPUT_PATH
  LIBRARY_OUTPUT_PATH)

include ( CheckIncludeFile )
include ( CheckFunctionExists )
include ( CheckTypeSize )

# Compiler defines... this should really be in a config file.
set( MOAB_DEFINES "" )
set( MOAB_LIBS )
set( MOAB_INSTALL_TARGETS )

################################################################################
# Options that the user controls
################################################################################
option ( BUILD_SHARED_LIBS   "Should shared or static libraries be created?"   ON  )
option ( ENABLE_CGM        "Should build with CGM support?"                  OFF )
option ( ENABLE_CGNS       "Should build with CGNS support?"                 OFF )
option ( ENABLE_MPI        "Should MOAB be compiled with MPI support?"       OFF )
option ( ENABLE_HDF5       "Include HDF I/O interfaces in the build?"                   OFF )
option ( ENABLE_ZLIB       "Include Zlib compression libraries (optionally used in HDF5)?"  OFF )
option ( ENABLE_SZIP       "Should build with szip support?"                 OFF )
option ( ENABLE_NETCDF     "Include NetCDF (ExodusII) interfaces in the build?" OFF )
option ( ENABLE_PNETCDF    "Include PNetCDF interfaces in the build (Requires NetCDF) ?" OFF )
option ( ENABLE_METIS      "Include Metis support for partitioning algorithms?" OFF )
option ( ENABLE_PARMETIS   "Include ParMetis support for partitioning algorithms?" OFF )
option ( ENABLE_ZOLTAN     "Include Zoltan support for partitioning algorithms?" OFF )
option ( ENABLE_VTK        "Include VTK I/O interfaces in the build?" OFF )
option ( ENABLE_TESTING "Enable Testing"                                  ON  )
option ( MOAB_FORCE_64_BIT_HANDLES "Force MBEntityHandle to be 64 bits (uint64_t)" OFF )
option ( MOAB_FORCE_32_BIT_HANDLES "Force MBEntityHandle to be 32 bits (uint32_t)" OFF )
option ( ENABLE_PYMOAB  "Build python interface to MOAB?" OFF )
option ( ENABLE_IMESH        "Should build IMESH?"       ON )
option ( ENABLE_FBIGEOM      "Should build FBIGEOM?"     OFF )
option ( ENABLE_IREL         "Should build IREL?"        OFF )
option ( ENABLE_MESQUITE     "Should build Mesquite mesh optimization interfaces?" OFF )
option ( ENABLE_BLASLAPACK   "Should require BLAS/LAPACK?" ON )
option ( ENABLE_FORTRAN      "Should enable FORTRAN?" ON )

mark_as_advanced(
  MOAB_FORCE_64_BIT_HANDLES
  MOAB_FORCE_32_BIT_HANDLES
)

if ( ENABLE_IMESH OR ENABLE_FBIGEOM )
  include (${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)
endif ()

################################################################################
# Check for system include files
################################################################################
check_include_file( inttypes.h   MOAB_HAVE_INTTYPES_H )
check_include_file( stdint.h     MOAB_HAVE_STDINT_H )
check_include_file( stddef.h     MOAB_HAVE_STDDEF_H )
check_include_file( stdlib.h     MOAB_HAVE_STDLIB_H )
check_include_file( sys/types.h  MOAB_HAVE_SYS_TYPES_H )
set( HAVE_SYS_TYPES_H ${MOAB_HAVE_SYS_TYPES_H} )
set( HAVE_STDDEF_H    ${MOAB_HAVE_STDDEF_H} )
set( HAVE_STDINT_H    ${MOAB_HAVE_STDINT_H} )
set( HAVE_INTTYPES_H    ${MOAB_HAVE_INTTYPES_H} )
set( HAVE_STDLIB_H    ${MOAB_HAVE_STDLIB_H} )
check_include_file( memory.h     MOAB_HAVE_MEMORY_H )

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)
################################################################################
# Integer size Related Settings
################################################################################
if ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )
  message( FATAL_ERROR
      "You may not turn both MOAB_FORCE_64_BIT_HANDLES and MOAB_FORCE_32_BIT_HANDLES on. Turn one off to continue."
    )
endif ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )

if ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )
  if ( MOAB_HAVE_INTTYPES_H )
    set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};inttypes.h" )
  endif ( MOAB_HAVE_INTTYPES_H )
  if ( MOAB_HAVE_STDLIB_H )
    set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};stdlib.h" )
    #set ( CHECK_TYPE_SIZE_PREMAIN "${CHECK_TYPE_SIZE_PREMAIN}\n#include <stdlib.h>\n" )
  endif ( MOAB_HAVE_STDLIB_H )
  check_type_size(  size_t       MOAB_HAVE_SIZE_T )
  check_type_size(  ptrdiff_t    MOAB_HAVE_PTRDIFF_T )
  set ( HAVE_SIZE_T ${MOAB_HAVE_SIZE_T} )
  set ( HAVE_PTRDIFF_T ${MOAB_HAVE_PTRDIFF_T} )
endif ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )

# Check of presence of finite or isfinite macros
include (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(std::isfinite MOAB_HAVE_STDISFINITE)
CHECK_FUNCTION_EXISTS(isfinite MOAB_HAVE_ISFINITE)
CHECK_FUNCTION_EXISTS(finite MOAB_HAVE_FINITE)

################################################################################
# Find packages
################################################################################

if ( BUILD_SHARED_LIBS EQUAL OFF )
#Static start
set_target_properties(MOAB PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(MOAB PROPERTIES LINK_SEARCH_END_STATIC 1)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif ( BUILD_SHARED_LIBS EQUAL OFF )

if ( ENABLE_BLASLAPACK  )
  # Find BLAS
  FIND_PACKAGE(BLAS REQUIRED)
  set (MOAB_HAVE_BLAS OFF CACHE INTERNAL "Configured MOAB with BLAS library support." )
  if (BLAS_FOUND)
    set(MOAB_HAVE_BLAS ON)
  endif(BLAS_FOUND)

  # Find LAPACK
  FIND_PACKAGE(LAPACK REQUIRED)
  set (MOAB_HAVE_LAPACK OFF CACHE INTERNAL "Configured MOAB with LAPACK library support." )
  if (LAPACK_FOUND)
    set(MOAB_HAVE_LAPACK ON)
  endif(LAPACK_FOUND)

  SET (ENABLE_FORTRAN ON)

  # Find Eigen3 is optional
  FIND_PACKAGE(Eigen3 OPTIONAL_COMPONENTS)
else (ENABLE_BLASLAPACK)
  # Find Eigen3 is required
  FIND_PACKAGE(Eigen3 REQUIRED)
endif (ENABLE_BLASLAPACK)

set (MOAB_HAVE_EIGEN OFF CACHE INTERNAL "Configured MOAB with Eigen3 library support." )
if ( EIGEN3_FOUND )
  INCLUDE_DIRECTORIES ( "${EIGEN3_DIR}" )
  SET(MOAB_DEFINES "${MOAB_DEFINES} ${EIGEN3_INCLUDES}")
  SET(MOAB_HAVE_EIGEN ON)
endif (EIGEN3_FOUND)

find_package( verdict REQUIRED )

# check for MPI package
set (MOAB_HAVE_MPI OFF CACHE INTERNAL "Found necessary MPI components. Configure MOAB with it." )
if ( ENABLE_MPI )
  set (MPI_C_COMPILER ${CC} )
  set (MPI_CXX_COMPILER ${CXX} )
  set (MPI_Fortran_COMPILER ${FC} )
  find_package( MPI REQUIRED )
  # CMake FindMPI script is sorely lacking:
  if ( MPI_LIBRARY AND MPI_INCLUDE_PATH )
    set( MPI_FOUND ON )
  endif ( MPI_LIBRARY AND MPI_INCLUDE_PATH )

  if ( MPI_FOUND )
    set ( MOAB_DEFINES "${MOAB_DEFINES} -DUSE_MPI" )
    set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    set( MOAB_HAVE_MPI ON )
  endif ( MPI_FOUND )
endif ( ENABLE_MPI )
include (config/CheckCompilerFlags.cmake)

set (MOAB_HAVE_ZLIB OFF CACHE INTERNAL "Found necessary Zlib components. Configure MOAB with it." )
if ( ENABLE_ZLIB )
  find_package( ZLIB REQUIRED )
  set (MOAB_HAVE_ZLIB ON)
endif (ENABLE_ZLIB)

set (MOAB_HAVE_SZIP OFF CACHE INTERNAL "Found necessary Zlib components. Configure MOAB with it." )
if ( ENABLE_SZIP )
  find_package( SZIP REQUIRED )
  set (MOAB_HAVE_SZIP ON)
endif (ENABLE_SZIP)

find_package( CURL )

set (MOAB_HAVE_HDF5 OFF CACHE INTERNAL "Found necessary HDF5 components. Configure MOAB with it." )
set (MOAB_HAVE_HDF5_PARALLEL OFF CACHE INTERNAL "Found necessary parallel HDF5 components. Configure MOAB with it." )
if ( ENABLE_HDF5 )
  if (MOAB_HAVE_MPI)
    set(HDF5_PREFER_PARALLEL TRUE)
  endif()
  find_package( HDF5_MOAB REQUIRED )
  set (MOAB_HAVE_HDF5 ON)
  # if ( HDF5_IS_PARALLEL AND MOAB_HAVE_MPI )
  if ( MOAB_HAVE_MPI )
    set (MOAB_HAVE_HDF5_PARALLEL ON)
  endif ( )
  if ( CURL_FOUND )
    set( HDF5_LIBRARIES ${HDF5_LIBRARIES} ${CURL_LIBRARIES} )
  endif ( CURL_FOUND )
  set( MOAB_LIBS ${HDF5_LIBRARIES} ${MOAB_LIBS} )
  include_directories( ${HDF5_INCLUDE_DIR} src/io/mhdf/include )
endif ( ENABLE_HDF5 )

set (MOAB_HAVE_NETCDF OFF CACHE INTERNAL "Found necessary NetCDF components. Configure MOAB with it." )
set (MOAB_HAVE_PNETCDF OFF CACHE INTERNAL "Found necessary PNetCDF components. Configure MOAB with it." )
if ( ENABLE_NETCDF )
  find_package( NetCDF REQUIRED )
  # set( MOAB_DEFINES "-DNETCDF_FILE ${MOAB_DEFINES}" )
  set (MOAB_HAVE_NETCDF ON)
  set (MOAB_HAVE_NETCDF_H 1)
  include_directories( ${NETCDF_INCLUDES} )
  set( MOAB_LIBS ${NETCDF_LIBRARIES} ${MOAB_LIBS} )
  #
  if ( PNETCDF_FOUND )
    # find_package( PNetCDF REQUIRED )
    set (MOAB_HAVE_PNETCDF ON)
    include_directories( ${PNETCDF_INCLUDES} )
    set( MOAB_LIBS ${PNETCDF_LIBRARIES} ${MOAB_LIBS} )
  endif ( PNETCDF_FOUND )
endif ( ENABLE_NETCDF )

set (MOAB_HAVE_METIS OFF CACHE INTERNAL "Found necessary Metis components. Configure MOAB with it." )
if ( ENABLE_METIS )
  find_package( Metis REQUIRED )
  set (MOAB_HAVE_METIS ON)
endif (ENABLE_METIS )

set (MOAB_HAVE_PARMETIS OFF CACHE INTERNAL "Found necessary ParMetis components. Configure MOAB with it." )
if ( ENABLE_PARMETIS AND ENABLE_METIS AND METIS_FOUND )
  find_package( ParMetis REQUIRED )
  set (MOAB_HAVE_PARMETIS ON)
endif ( ENABLE_PARMETIS AND ENABLE_METIS AND METIS_FOUND )

set (MOAB_HAVE_ZOLTAN OFF CACHE INTERNAL "Found necessary Zoltan components. Configure MOAB with it." )
if ( ENABLE_ZOLTAN )
  find_package( Zoltan REQUIRED )
  set (MOAB_HAVE_ZOLTAN ON)
endif (ENABLE_ZOLTAN )

set (MOAB_HAVE_CGM OFF CACHE INTERNAL "Found necessary CGM components. Configure MOAB with it." )
if ( ENABLE_CGM )
   find_package( CGM REQUIRED )
   set( MOAB_DEFINES "${CGM_DEFINES} -DCGM ${MOAB_DEFINES}" )
   set (MOAB_HAVE_CGM ON)
   set (ENABLE_IREL ON)
   set (ENABLE_FBIGEOM ON)
endif (ENABLE_CGM)

set (MOAB_HAVE_CGNS OFF CACHE INTERNAL "Found necessary CGNS components. Configure MOAB with it." )
if (ENABLE_CGNS)
  set( MOABIO_LIBS ${MOABIO_LIBS} ${CGNS_LIBRARIES} )
  set (MOAB_HAVE_CGNS ON)
endif(ENABLE_CGNS)

# Strip our custom flags
string(STRIP "${CMAKE_EXE_LINKER_FLAGS} " CMAKE_EXE_LINKER_FLAGS)
string(STRIP "${MOAB_CXX_FLAGS} " MOAB_CXX_FLAGS)
string(STRIP "${MOAB_DEFINES} " MOAB_DEFINES)

################################################################################
# Add Directories
################################################################################
add_subdirectory( src )
add_subdirectory( itaps )
add_subdirectory( tools )
if (ENABLE_PYMOAB)
  add_subdirectory( pymoab )
endif()
################################################################################
# Testing Related Settings
################################################################################
#turn on ctest if we want testing
if ( ENABLE_TESTING )
  enable_testing()
  add_subdirectory( test )
endif()

###############################################################################
#
###############################################################################
export(
  TARGETS ${MOAB_INSTALL_TARGETS}
  FILE "${PROJECT_BINARY_DIR}/MOABTargets.cmake"
  )
install(
  EXPORT MOABTargets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB/
  )

# Create some custom command

# For consistency with autoconf, create a "make check" equivalent to "ctest"
ADD_CUSTOM_TARGET( check COMMAND ctest DEPENDS )

# For consistency with autoconf, create "make dist" and "make distcheck" targets
INCLUDE( CMakeDistCheck )
DEFINE_DISTCHECK ( )

################################################################################
# Generate the MOABConfig.cmake file
################################################################################
set(CXX ${CMAKE_CXX_COMPILER})
set(CC ${CMAKE_C_COMPILER})
set(F77 ${CMAKE_Fortran_COMPILER} )
set(FC ${CMAKE_Fortran_COMPILER})
set(CPPFLAGS "${CMAKE_CXX_FLAGS} ${MOAB_DEFINES}")
IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG})
  set(CFLAGS ${CMAKE_C_FLAGS_DEBUG})
  set(FFLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
ELSE (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE})
  set(CFLAGS ${CMAKE_C_FLAGS_RELEASE})
  set(FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
ENDIF (CMAKE_BUILD_TYPE MATCHES "Debug")

# CMake does not define uninstall target. Look at install manifest and remove
# all the files recursively.
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
    )

if (ENABLE_FORTRAN)
  enable_language(Fortran)

  # Include CMake's Fortran mangling evaluation macros
  include(FortranCInterface)

  if ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
    set(FC_CASE "NAME")
  else ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
    set(FC_CASE "name")
  endif ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
  if (NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
    set(FortranCInterface_GLOBAL_SUFFIX "## ${FortranCInterface_GLOBAL_SUFFIX}")
  endif (NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
  if (NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")
    set(FortranCInterface_GLOBAL__SUFFIX "## ${FortranCInterface_GLOBAL__SUFFIX}")
  endif (NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")

  set(MOAB_F77_FUNC  "${FortranCInterface_GLOBAL_PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL_SUFFIX}")
  set(MOAB_F77_FUNC_ "${FortranCInterface_GLOBAL__PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL__SUFFIX}")
  set(MOAB_FC_FUNC   "${FortranCInterface_GLOBAL_PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL_SUFFIX}")
  set(MOAB_FC_FUNC_  "${FortranCInterface_GLOBAL__PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL__SUFFIX}")
endif (ENABLE_FORTRAN)

# Generate the CMake input header file and then pass it through configuration
## execute_process( COMMAND cat MOABConfig.h.in COMMAND sed -e "s/#undef /#cmakedefine MOAB_/g"
#execute_process( COMMAND cat MOABConfig.h.in COMMAND sed -r "s/#undef ([A-Z].*)/#cmakedefine MOAB_\\1 @MOAB_\\1@/"
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/config
#  OUTPUT_FILE MOABConfig.cmake2.h.in
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#  )
CONFIGURE_FILE(config/MOABConfig.cmake.h.in "${PROJECT_BINARY_DIR}/src/moab/MOABConfig.h" @ONLY)

CONFIGURE_FILE(config/MOABConfig.cmake.cmake "${PROJECT_BINARY_DIR}/MOABConfig.cmake" @ONLY)
CONFIGURE_FILE(moab.make.in "${PROJECT_BINARY_DIR}/lib/moab.make" @ONLY)
if ( MOAB_HAVE_FBIGEOM )
  CONFIGURE_FILE(itaps/fbigeom/FBiGeom-Defs.inc.in "${PROJECT_BINARY_DIR}/lib/FBiGeom-Defs.inc" @ONLY)
  install(
    FILES "${PROJECT_BINARY_DIR}/lib/FBiGeom-Defs.inc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( MOAB_HAVE_FBIGEOM )

if ( MOAB_HAVE_IMESH )
  CONFIGURE_FILE(itaps/imesh/iMesh-Defs.inc.in "${PROJECT_BINARY_DIR}/lib/iMesh-Defs.inc" @ONLY)
  INSTALL(
    FILES "${PROJECT_BINARY_DIR}/lib/iMesh-Defs.inc"
    "${PROJECT_BINARY_DIR}/lib/moab.make"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( MOAB_HAVE_IMESH )
INSTALL(
  FILES "${PROJECT_BINARY_DIR}/MOABConfig.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB/
  )
INSTALL(
  FILES "${PROJECT_SOURCE_DIR}/config/ResolveCompilerPaths.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB
  )

ADD_CUSTOM_TARGET(
  uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  )

# All done.
