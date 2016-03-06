#######################################################################################
# Check values of NC_MAX_DIMS and NC_MAX_VARS for ExodusII compatability.
# 
# Arguments are: 
#  1) required value for NC_MAX_DIMS, may be emtpy
#  2) required value for NC_MAX_VARS, may be emtpy
#  3) name of header in which to check for values
#  4) variable in which to store result (set to yes if
#        limits are sufficient, and on otherwise.)
#######################################################################################
AC_DEFUN([FATHOM_CHECK_NETCDF_LIMITS],[
  MIN_NC_MAX_DIMS="$1"
  MIN_NC_MAX_VARS="$2"
  $4=yes

  if test "x" != "x$MIN_NC_MAX_DIMS"; then
    AC_MSG_CHECKING([if NC_MAX_DIMS is at least ${MIN_NC_MAX_DIMS}])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([#include <$3>],
                       [[int arr[1 + (int)(NC_MAX_DIMS) - (int)(${MIN_NC_MAX_DIMS})];]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no]); $4=no])
  fi
  if test "x" != "x$MIN_NC_MAX_VARS"; then
    AC_MSG_CHECKING([if NC_MAX_VARS is at least ${MIN_NC_MAX_VARS}])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([#include <netcdf.h>],
                       [[int arr[1 + (int)(NC_MAX_VARS) - (int)(${MIN_NC_MAX_VARS})];]])],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no]); $4=no])
  fi
])

#######################################################################################
# Check for NetCDF library
# Sets HAVE_NETCDF to 'yes' or 'no'
# If HAVE_NETCDF == yes, then exports:
#   NETCDF_CPPFLAGS
#   NETCDF_LDFLAGS
#   NETCDF_LIBS
#   NETCDF_SUFFICIENT_DIMS_VARS
#
# This macro has two optional arguments:  a minimum value for
# NC_MAX_DIMS and a minimum value for NC_MAX_VARS.  If either or
# both of these are specified, NETCDF_SUFFICIENT_DIMS_VARS will
# be set to yes if the NetCDF library is built with limits that
# are at least the passed miminums.  It will be set to no if 
# either limit is less than the passed minimum.
#######################################################################################
AC_DEFUN([FATHOM_CHECK_NETCDF],[

  # Supported NetCDF versions: 4.3.3, 4.3.2, 4.2.1
  # Arguments: 1) Default Version Number, 2) Download by default ?
  AUSCM_CONFIGURE_DOWNLOAD_NETCDF([4.3.2],[no])

AC_MSG_CHECKING([if NetCDF support is enabled])
AC_ARG_WITH(netcdf, 
[AS_HELP_STRING([--with-netcdf@<:@=DIR@:>@], [Specify NetCDF library to use for ExodusII file format])
AS_HELP_STRING([--without-netcdf], [Disable support for ExodusII file format])],
[NETCDF_DIR=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-netcdf=\"${withval}\""
]
, [NETCDF_DIR=$NETCDF_DIR])
if (test "x" != "x$NETCDF_DIR" && test "xno" != "x$NETCDF_DIR"); then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

# if NetCDF support is not disabled
HAVE_NETCDF=no
if (test "x" != "x$NETCDF_DIR" && test "xno" != "x$NETCDF_DIR"); then

    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$NETCDF_DIR" && test "x" != "x$NETCDF_DIR"; then
    if test -d "${NETCDF_DIR}/lib"; then
      NETCDF_LDFLAGS="-L${NETCDF_DIR}/lib"
    elif test -d "${NETCDF_DIR}"; then
      NETCDF_LDFLAGS="-L${NETCDF_DIR}"
    else
      AC_MSG_ERROR("$NETCDF_DIR is not a directory.")
    fi
    if test -d "${NETCDF_DIR}/include"; then
      NETCDF_CPPFLAGS="-I${NETCDF_DIR}/include"
    elif test -d "${NETCDF_DIR}/inc"; then
      NETCDF_CPPFLAGS="-I${NETCDF_DIR}/inc"
    else
      NETCDF_CPPFLAGS="-I${NETCDF_DIR}"
    fi
  fi
 
  AC_PATH_PROG([NC_CONFIG], [nc-config], [no], [${NETCDF_DIR}/bin])
  if test "$NC_CONFIG" != "no" ; then
    NETCDF_CPPFLAGS="`$NC_CONFIG --cflags`"
    NETCDF_LDFLAGS="`$NC_CONFIG --libs`"
    NETCDF_VERSION="`$NC_CONFIG --version`"
    AC_SUBST([NETCDF_CPPFLAGS])
    AC_SUBST([NETCDF_LDFLAGS])
    AC_SUBST([NETCDF_VERSION])
    HAVE_NETCDF=yes
  else
    AC_MSG_WARN("NetCDF configuration utility not found")
  fi
   
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$NETCDF_CPPFLAGS $CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$NETCDF_LDFLAGS $PNETCDF_LDFLAGS $LDFLAGS"
  
   # Check for C library
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS( [netcdf.h], 
                    [FATHOM_CHECK_NETCDF_LIMITS([$1],[$2],[netcdf.h],[HAVE_NETCDF=yes; NETCDF_SUFFICIENT_DIM_VARS])], 
                    [AC_MSG_WARN([[NetCDF header not found. Disabling NetCDF configuration.]]); HAVE_NETCDF=no; NETCDF_DIR="";] )

  # Check if netcdf is usable by itself
  if test "x$HAVE_NETCDF" != "xno"; then
    AC_CHECK_LIB( [netcdf], [nc_create], [NETCDF_LIBS="-lnetcdf $PNETCDF_LIBS"; HAVE_NETCDF=yes;], [
      # Check if netcdf is usable with HDF5
      unset ac_cv_lib_netcdf
      unset ac_cv_lib_netcdf_nc_create
      # If we haven't already looked for HDF5 libraries, again now incase
      # they're in the NetCDF lib directory.
      FATHOM_DETECT_HDF5_LIBS
      LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
      AC_CHECK_LIB( [netcdf], [nc_create], [NETCDF_LIBS="-lnetcdf $PNETCDF_LIBS -lhdf5_hl $HDF5_LIBS"; HAVE_NETCDF=yes;], [
        # Try one more time with HDF5 and libcurl
        unset ac_cv_lib_netcdf
        unset ac_cv_lib_netcdf_nc_create
        AC_CHECK_LIB( [netcdf], [nc_create], [NETCDF_LIBS="-lnetcdf $PNETCDF_LIBS -lhdf5_hl $HDF5_LIBS -lcurl"; HAVE_NETCDF=yes;],
          [HAVE_NETCDF=no], [$PNETCDF_LIBS -lhdf5_hl $HDF5_LIBS -lcurl ] )],
        [$PNETCDF_LIBS -lhdf5_hl $HDF5_LIBS] )],
      )
  fi

  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  AC_LANG_POP([C])

  if test "x$HAVE_NETCDF" = "xno"; then
    if test "x$NETCDF_DIR" != "x"; then 
      AC_MSG_ERROR("NetCDF not found or not working")
    else
      AC_MSG_WARN("NetCDF support disabled")
    fi
    NETCDF_CPPFLAGS=
    NETCDF_LDFLAGS=
  fi
fi

]) # FATHOM_HAVE_NETCDF

AC_DEFUN([FATHOM_CHECK_PNETCDF],[

AC_ARG_WITH(pnetcdf, 
[AS_HELP_STRING([--with-pnetcdf@<:@=DIR@:>@], [Specify PNetCDF library to use])
AS_HELP_STRING([--without-pnetcdf], [Disable support for PNetCDF-based file formats])],
[PNETCDF_DIR=$withval
DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-pnetcdf=\"${withval}\""
]
, [PNETCDF_DIR=])

  # PNETCDF requires MPI too
if test "xyes" != "x$enablempi"; then
  if test "x" == "x$PNETCDF_DIR"; then
    PNETCDF_DIR=no
  elif test "xno" != "xPNETCDF_DIR"; then
    AC_MSG_ERROR([PNetCDF requires --with-mpi])
  fi
  HAVE_PNETCDF=no
fi

AC_MSG_CHECKING([if PNetCDF support is enabled])
if (test "xno" != "x$PNETCDF_DIR" && test "x$PNETCDF_DIR" != "x"); then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

 # if Pnetcdf support is not disabled
HAVE_PNETCDF=no
if (test "xno" != "x$PNETCDF_DIR" && test "x$PNETCDF_DIR" != "x"); then
  HAVE_PNETCDF=yes
  
    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$PNETCDF_DIR" && test "x" != "x$PNETCDF_DIR"; then
    if test -d "${PNETCDF_DIR}/lib"; then
      PNETCDF_LDFLAGS="-L${PNETCDF_DIR}/lib"
    elif test -d "${PNETCDF_DIR}"; then
      PNETCDF_LDFLAGS="-L${PNETCDF_DIR}"
    else
      AC_MSG_ERROR("$PNETCDF_DIR is not a directory.")
    fi
    if test -d "${PNETCDF_DIR}/include"; then
      PNETCDF_CPPFLAGS="-I${PNETCDF_DIR}/include"
    elif test -d "${PNETCDF_DIR}/inc"; then
      PNETCDF_CPPFLAGS="-I${PNETCDF_DIR}/inc"
    else
      PNETCDF_CPPFLAGS="-I${PNETCDF_DIR}"
    fi
  fi
  
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$PNETCDF_CPPFLAGS $CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$PNETCDF_LDFLAGS $LDFLAGS"
  
   # Check for C library
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS( [pnetcdf.h], 
                    [FATHOM_CHECK_NETCDF_LIMITS([$1],[$2],[pnetcdf.h],[PNETCDF_SUFFICIENT_DIM_VARS])], 
                    [HAVE_PNETCDF=no] )

  AC_CHECK_LIB( [pnetcdf], [ncmpi_create], [PNETCDF_LIBS="-lpnetcdf"], [HAVE_PNETCDF=no] )

  
  AC_LANG_POP([C])
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"

  if test "x$HAVE_PNETCDF" = "xno"; then
    if test "x$PNETCDF_DIR" != "x"; then 
      AC_MSG_ERROR("PNetCDF not found or not working")
    fi
    PNETCDF_CPPFLAGS=
    PNETCDF_LDFLAGS=
  fi
fi

]) # FATHOM_HAVE_PNETCDF
