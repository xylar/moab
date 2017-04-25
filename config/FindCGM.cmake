# FindCGM.cmake
#
# If you set the CGM_CFG CMake variable to point to a file named "cgm.make"
# (produced by CGM as part of any build or install), then the script will
# locate CGM assets for your package to use.
#
# This script defines the following CMake variables:
#   CGM_FOUND         defined when CGM is located, false otherwise
#   CGM_INCLUDE_DIRS  directories containing CGM headers
#   CGM_DEFINES       preprocessor definitions you should add to source files
#   CGM_LIBRARIES     paths to CGM library and its dependencies
set (CGM_DIR "" CACHE PATH "Path to search for CGM header and library files")
find_file(CGM_CMAKE_CFG CGMConfig.cmake
    HINTS ${CGM_DIR} ${CGM_DIR}/lib/cmake/CGM
)
MARK_AS_ADVANCED(CGM_CMAKE_CFG)

message("CGM config path - ${CGM_CMAKE_CFG}")
if (CGM_CMAKE_CFG)
include (${CGM_CMAKE_CFG})
endif (CGM_CMAKE_CFG)

if (CGM_FOUND)

# Include dir
find_path(IGEOM_INCLUDE_DIR
    NAMES iGeom.h
    PATHS ${CGM_INCLUDES}
)

find_path(CGM_LIB_PATH
    NAMES cgm.make
    PATHS ${CGM_DIR} ${CGM_DIR}/lib ${CGM_DIR}/../..
)
MARK_AS_ADVANCED(IGEOM_INCLUDE_DIR CGM_LIB_PATH)
#get_filename_component(OCE_LIBRARY_DIR "${OCE_DIR}/OCEConfig.cmake" PATH)

SET(MOAB_HAVE_CGM_FACET ${CGM_HAS_FACET})
SET(MOAB_HAVE_CGM_OCC ${CGM_HAS_OCC})
SET(MOAB_CGM_GEOM_ENGINE "${PRIMARY_GEOMETRY_ENGINE}")

SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CGM_LDFLAGS}")

# Split the version correctly
string(REPLACE "." ";" VERSION_LIST "${CGM_VERSION}")
list(GET VERSION_LIST 0 CGM_MAJOR_VERSION)
list(GET VERSION_LIST 1 CGM_MINOR_VERSION)

# Output details about CGM configuration
message("Found CGM: ${CGM_DIR}")
message("---   CGM configuration ::")
message("        Primary Geometry Engine : ${PRIMARY_GEOMETRY_ENGINE}")
message("        Include Directory       : ${IGEOM_INCLUDE_DIR}")
message("        Library Directory       : ${CGM_LIB_PATH}")

endif(CGM_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CGM
  REQUIRED_VARS CGM_INCLUDES CGM_LDFLAGS CGM_LIBRARIES IGEOM_INCLUDE_DIR
  VERSION_VAR CGM_VERSION
)
