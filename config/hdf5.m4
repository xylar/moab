AC_DEFUN([FATHOM_HDF5_LIBS_HELPER],[
if (test "x$HAVE_LIB_HDF5" != "xyes"); then
   unset "ac_cv_lib_${HDF5_LIBNAME}_H5Fopen"
   unset "ac_cv_lib_${HDF5_LIBNAME}___H5Fopen"
   AC_CHECK_LIB( [${HDF5_LIBNAME}], [H5Fopen], [HAVE_LIB_HDF5=yes; HDF5_LIBS="$HDF5_LIBS -l${HDF5_LIBNAME} $1"], [], [$1] )
fi
])

dnl ---------------------------------------------------------------------------
dnl FATHOM_HDF5_LIBS_HELPER_FORTRAN
dnl   Inserts the correct fortran libraries into HDF5 libraries
dnl   Arguments:
dnl ---------------------------------------------------------------------------
AC_DEFUN([FATHOM_HDF5_LIBS_HELPER_FORTRAN],[
  # Make Fortran link line by inserting Fortran libraries
  for arg in $HDF5_LIBS
  do
    case "$arg" in
      -lhdf5_hl) AC_HAVE_LIBRARY([hdf5hl_fortran], [HDF5_FLIBS="$HDF5_FLIBS -lhdf5hl_fortran"], [], [$HDF5_LIBS $HDF5_FLIBS])
        ;; 
      -lhdf5)    AC_HAVE_LIBRARY([hdf5_fortran], [HDF5_FLIBS="$HDF5_FLIBS -lhdf5_fortran"], [], [$HDF5_LIBS $HDF5_FLIBS])
        ;; 
    esac
  done
  HDF5_LIBS="$HDF5_LIBS $HDF5_FLIBS"
])


dnl ---------------------------------------------------------------------------
dnl FATHOM_HDF5_LIBS_HELPER_CXX
dnl   Inserts the correct C++ libraries into HDF5 libraries
dnl   Arguments:
dnl ---------------------------------------------------------------------------
AC_DEFUN([FATHOM_HDF5_LIBS_HELPER_CXX],[
  # Make CXX link line by inserting C++ libraries
  oldLIBS=$LIBS
  LIBS=""
  for arg in $HDF5_LIBS
  do
    case "$arg" in #(
      -lhdf5_hl) AC_HAVE_LIBRARY([hdf5hl_cpp], [HDF5_CXXLIBS="$HDF5_CXXLIBS -lhdf5hl_cpp"], [], [$HDF5_LIBS $HDF5_CXXLIBS])
        ;; #(
      -lhdf5)    AC_HAVE_LIBRARY([hdf5_cpp], [HDF5_CXXLIBS="$HDF5_CXXLIBS -lhdf5_cpp"], [], [$HDF5_LIBS $HDF5_CXXLIBS])
        ;; #(
    esac
  done
  HDF5_LIBS="$HDF5_LIBS $HDF5_CXXLIBS"
  LIBS=$oldLIBS
])


dnl ---------------------------------------------------------------------------
dnl FATHOM_HDF5_LIBS_HELPER_HL
dnl   Inserts the correct high level libraries into HDF5 libraries
dnl   Arguments:
dnl ---------------------------------------------------------------------------
AC_DEFUN([FATHOM_HDF5_LIBS_HELPER_HL],[
if (test $HAVE_LIB_HDF5 = yes && test $HAVE_LIB_HDF5HL = no); then
  # Look for HDF5's high level library
  AC_HAVE_LIBRARY([hdf5_hl], [HAVE_LIB_HDF5HL=yes; HDF5_LIBS="-lhdf5_hl $HDF5_LIBS"], [], [$HDF5_LIBS])
  if (test $HAVE_LIB_HDF5HL=yes); then
    HDF5_FLIBS=""
    FATHOM_HDF5_LIBS_HELPER_FORTRAN
    HDF5_CXXLIBS=""
    FATHOM_HDF5_LIBS_HELPER_CXX
  fi
fi
])

#######################################################################################
# Helper function for FATHOM_CHECK_HDF5 and FATHOM_CHECK_NETCDF
# If HAVE_LIB_HDF5 == yes, then does nothing.
# Otherwise sets HAVE_LIB_HDF5 to yes or no depending on whether or not
# HDF5 libraries are detected and sets HDF5_LIBS
# Respects caller's LDFLAGS, CPPFLAGS, and LIBS.  Caller should set appropriately.
# Respections optional HDF5_LIBNAME uses to specify alternate library name.  If
# not specified, will be set to -lhdf5
#######################################################################################
AC_DEFUN([FATHOM_DETECT_HDF5_LIBS],[

 # if we've already done this check, then don't do it again
if (test "xyes" != "x$HAVE_LIB_HDF5"); then
  test "x" != "x$HDF5_LIBNAME" || HDF5_LIBNAME=hdf5
  
  HAVE_LIB_HDF5=no
  HAVE_LIB_HDF5HL=no
  FATHOM_HDF5_LIBS_HELPER([$LIBS])
  FATHOM_HDF5_LIBS_HELPER_HL
fi
])

#######################################################################################
# Check for HDF5 library and related stuff
# Sets enablehdf5 to 'yes' or 'no'
# If enablehdf5 == yes, then sets:
#  HDF5_CPPFLAGS
#  HDF5_LDFLAGS
#  HDF5_LIBS
#######################################################################################
AC_DEFUN([FATHOM_CHECK_HDF5],[
  
  # CLI option for linking zlib
AC_ARG_WITH(zlib,
  [AS_HELP_STRING([--with-zlib=DIR],[HDF5 requires zlib, and zlib can be found at...])],
  [if (test "x$withval" != "x" && test "x$withval" != "xno"); then 
    ZLIB_DIR=$withval
    DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-zlib=\"${withval}\""
  fi
  ],[ZLIB_DIR=])
case "x$ZLIB_DIR" in
  xyes|xno|x)
    ;;
  *)
    if ! test -d  ${ZLIB_DIR}/lib; then
      AC_MSG_ERROR([Not a directory: ${ZLIB_DIR}/lib])
    fi
    HDF5_LDFLAGS="$HDF5_LDFLAGS -L${ZLIB_DIR}/lib"
    ;;
esac
enablezlib=no
if test "x$ZLIB_DIR" != "xno"; then
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
  AC_CHECK_LIB([z],[deflate],[enablezlib=yes; HDF5_LIBS="$HDF5_LIBS -lz"],
    [if test "x$ZLIB_DIR" != "x"; then AC_MSG_ERROR([Could not find zlib]); fi])
  LDFLAGS="$old_LDFLAGS"
fi
AC_SUBST(enablezlib)
AC_SUBST(ZLIB_DIR)

  # CLI option for linking szip
AC_ARG_WITH(szip,
  [AS_HELP_STRING([--with-szip=DIR],[HDF5 requires szip, and szip an be found at...])],
  [if (test "x$withval" != "x" && test "x$withval" != "xno"); then
    SZIP_DIR=$withval
    DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-szip=\"${withval}\""
  fi],[SZIP_DIR=])
case "x$SZIP_DIR" in
  xyes|xno|x)
    ;;
  *)
    if ! test -d  ${SZIP_DIR}/lib; then
      AC_MSG_ERROR([Not a directory: ${SZIP_DIR}/lib])
    fi
    HDF5_LDFLAGS="$HDF5_LDFLAGS -L${SZIP_DIR}/lib"
    ;;
esac
enableszip=no
if test "x$SZIP_DIR" != "xno"; then
  old_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
  AC_CHECK_LIB([sz],[SZ_Decompress],[enableszip=yes; HDF5_LIBS="$HDF5_LIBS -lsz"],
    [if test "x$SZIP_DIR" != "x"; then AC_MSG_ERROR([Could not find libsz]); fi])
  LDFLAGS="$old_LDFLAGS"
fi
AC_SUBST(enableszip)
AC_SUBST(SZIP_DIR)

  # CLI option for extra HDF5 link flags
AC_ARG_WITH([hdf5-ldflags],[AS_HELP_STRING([--with-hdf5-ldflags=...],
 [Extra LDFLAGS required for HDF5 library (e.g. parallel IO lib)])],
 [HDF5_LDFLAGS_WITHVAL="$withval"; 
 DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5-ldflags=\"${withval}\""
],[HDF5_LDFLAGS_WITHVAL=])
case "x$HDF5_LDFLAGS_WITHVAL" in
  xno)
    AC_MSG_ERROR("Invalid argument: --without-hdf5-ldflags")
    ;;
  xyes)
    AC_MSG_ERROR("Nonsensical argument:  --with-hdf5-ldflags without any specified flags")
    ;;
  *)
    HDF5_LDFLAGS="$HDF5_LDFLAGS $HDF5_LDFLAGS_WITHVAL"
    ;;
esac


  # CLI option for HDF5
AC_MSG_CHECKING([if HDF5 support is enabled])
AC_ARG_WITH(hdf5, 
[AS_HELP_STRING([--with-hdf5@<:@=DIR@:>@], [Specify HDF5 library to use for native file format])
AS_HELP_STRING([--without-hdf5], [Disable support for native HDF5 file format])],
[if (test "x$withval" != "x" && test "x$withval" != "xno"); then
  HDF5_DIR=$withval
  DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hdf5=\"${withval}\""
 fi
], [HDF5_DIR=$HDF5_DIR])
if (test "x" != "x$HDF5_DIR" && test "xno" != "x$HDF5_DIR"); then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  # Reset the directory since we do not want to configure HDF5
  HDF5_DIR=""
fi

# Supported HDF5 versions: 1.8.10, 1.8.12, 1.8.14, 1.8.15
# Arguments: 1) Default Version Number, 2) Download by default ?
AUSCM_CONFIGURE_DOWNLOAD_HDF5([1.8.12],[no])

enablehdf5=no
if (test "x" != "x$HDF5_DIR" && test "xno" != "x$HDF5_DIR"); then
  enablehdf5=yes

    # if a path is specified, update LIBS and INCLUDES accordingly
  if test "xyes" != "x$HDF5_DIR" && test "x" != "x$HDF5_DIR"; then
    if test -d "${HDF5_DIR}/dll"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_DIR}/dll"
      HDF5_LIBNAME=hdf5dll
    elif test -d "${HDF5_DIR}/lib"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_DIR}/lib"
    elif test -d "${HDF5_DIR}"; then
      HDF5_LDFLAGS="$HDF5_LDFLAGS -L${HDF5_DIR}"
    else
      AC_MSG_ERROR("$HDF5_DIR is not a directory.")
    fi
    if test -d "${HDF5_DIR}/include"; then
      HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_DIR}/include"
      if test "x$GXX" = "xyes" && test "x$GCC" = "xyes"; then
        HDF5_CPPFLAGS="$HDF5_CPPFLAGS -isystem ${HDF5_DIR}/include"
      fi
    fi
    if test -d "${HDF5_DIR}/include/hdf5/include"; then
      HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_DIR}/include/hdf5/include"
    else
      if test -d "${HDF5_DIR}/include/hdf5"; then
        HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_DIR}/include/hdf5"
      else
        HDF5_CPPFLAGS="$HDF5_CPPFLAGS -I${HDF5_DIR}"
      fi
    fi
  fi
 
  # Check for IBM parallel IO library
  if test "x$enablempi" != "xno"; then
    AC_CHECK_LIB([gpfs],[gpfs_stat],[LIBS="-lgpfs $LIBS"])
  fi

  # Add flag to defines
  old_CPPFLAGS="$CPPFLAGS"
  old_LDFLAGS="$LDFLAGS"
  old_LIBS="$LIBS"
  CPPFLAGS="$HDF5_CPPFLAGS $CPPFLAGS"
  LDFLAGS="$HDF5_LDFLAGS $LDFLAGS"
  LIBS="$HDF5_LIBS $LIBS"
  
  # check for libraries and headers
  AC_CHECK_HEADERS( [hdf5.h], [], [enablehdf5=no] )

  HAVE_LIB_HDF5=no
  FATHOM_DETECT_HDF5_LIBS

  if (test "x$HAVE_LIB_HDF5" != "xyes"); then
    enablehdf5=no
  fi
  
  if (test "x$enablehdf5" != "xyes"); then
    if test "x" = "x$HDF5_DIR"; then
      AC_MSG_WARN([HDF5 library not found or not usable.])
    else
      AC_MSG_ERROR([HDF5 library not found or not usable.])
    fi
    HDF5_CPPFLAGS=
    HDF5_LDFLAGS=
    HDF5_LIBS=
  fi
  
  CPPFLAGS="$old_CPPFLAGS"
  LDFLAGS="$old_LDFLAGS"
  LIBS="$old_LIBS"
fi

if (test "xyes" != "x$enablehdf5"); then
  AC_MSG_WARN([Support for native HDF5 file format disabled])
else
  AC_DEFINE([HAVE_HDF5],[1],[Define if configured with HDF5 support.])
fi
AM_CONDITIONAL(HAVE_HDF5, [test "xno" != "x$enablehdf5"])
AC_SUBST(enablehdf5)
AC_SUBST(HDF5_DIR)
AC_SUBST(HDF5_CPPFLAGS)
MB_CPPFLAGS="$HDF5_CPPFLAGS $MB_CPPFLAGS"
EXPORT_LDFLAGS="$EXPORT_LDFLAGS $HDF5_LDFLAGS"
AC_SUBST(HDF5_LIBS)

WARN_PARALLEL_HDF5=no
WARN_PARALLEL_HDF5_NO_COMPLEX=no
enablehdf5parallel=no
if (test "xno" != "x$enablehdf5"); then
  if (test "xno" != "x$enablempi"); then
    old_LDFLAGS="$LDFLAGS"
    old_LIBS="$LIBS"
    LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
    LIBS="$HDF5_LIBS $LIBS"
    IS_HDF5_NOT_PARALLEL=""
    AC_PATH_PROGS([H5CC], [h5cc h5pcc], [], [$HDF5_DIR/bin])
    AC_MSG_CHECKING([for parallel HDF5 support])
    if (test "x$H5CC" != "x"); then dnl check if HDF5 was configured with parallel support
      IS_HDF5_NOT_PARALLEL="`$H5CC -showconfig | grep 'Parallel HDF5: no'`"
    fi
    if (test "x$IS_HDF5_NOT_PARALLEL" != "x"); then
      AC_MSG_RESULT(no)
      WARN_PARALLEL_HDF5=yes
      AC_MSG_WARN("libhdf5 library does not include parallel support.  Parallel HDF5 I/O disabled")
    else
      AC_MSG_RESULT(yes)
      enablehdf5parallel=yes
    fi

    LDFLAGS="$old_LDFLAGS"
    LIBS="$old_LIBS"
  fi
fi
AM_CONDITIONAL(HAVE_HDF5_PARALLEL, [test "xno" != "x$enablehdf5parallel"])
if test "xno" != "x$enablehdf5parallel"; then
  AC_DEFINE([HAVE_HDF5_PARALLEL],[1],[Define if configured with Parallel HDF5 support.])
  
  AC_MSG_CHECKING([for H5_MPI_COMPLEX_DERIVED_DATATYPE_WORKS])
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
  AC_PREPROC_IFELSE([AC_LANG_PROGRAM([#include <H5pubconf.h>],[
#ifndef H5_MPI_COMPLEX_DERIVED_DATATYPE_WORKS
  choke me
#endif])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no); WARN_PARALLEL_HDF5_NO_COMPLEX=yes])
  CPPFLAGS="$old_CPPFLAGS"
fi
AC_SUBST(enablehdf5parallel)

]) # FATHOM_CHECK_HDF5
