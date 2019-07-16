dnl ----------------------------------------------------------------
dnl Locate header files for the C++ linear algebra library Eigen.
dnl Eigen is a header-only template library. By default we check for the
dnl Eigen files in the --with-eigen-include=xxx argument provided to
dnl configure, or if those don't exist in the $EIGEN3_DIR/Eigen directory,
dnl or in /usr/include.  
dnl
dnl Note: Eigen is installed (by default) at the location
dnl /path/to/eigen/Eigen, i.e. with path ending in capital 'Eigen'.
dnl You should specify --with-eigen-include=/path/to/eigen
dnl during configure, or set your $EIGEN3_DIR environment variable
dnl to /path/to/eigen.
dnl
dnl - Modified lightly based on libMesh's version.
dnl ----------------------------------------------------------------

AC_DEFUN([FATHOM_CONFIGURE_EIGEN3], 
[
  dnl User-specific include path
  AC_ARG_WITH(eigen3,
              AC_HELP_STRING([--with-eigen3=PATH],[Specify the path for Eigen3 header files]),
              witheigeninc=$withval,
              witheigeninc=no)

  dnl Fall back on default paths to Eigen's include files
  if (test "x$witheigeninc" != "xno"); then
    EIGEN3_DIR="$witheigeninc"
  elif test "x$EIGEN3_DIR" != x -a -f $EIGEN3_DIR/Eigen/Eigen; then
    echo "Environment EIGEN3_DIR=$EIGEN3_DIR"
  elif test -f /usr/include/eigen3/Eigen/Eigen ; then
    EIGEN3_DIR="/usr/include/eigen3"
  elif test -f /usr/local/include/eigen3/Eigen/Eigen ; then
    EIGEN3_DIR="/usr/local/include/eigen3"
  else
    EIGEN3_DIR=""
  fi

  enableeigen=no;
  if (test "x$EIGEN3_DIR" != "x"); then
 
    dnl Check for existence of a header file in the specified location.  Note: here
    dnl we are checking for the header file "Eigen" in the Eigen directory.
    dnl AC_CHECK_FILE([$EIGEN3_DIR/Eigen], [eigenincFound="OK"], [eigenincFound="FAIL"])
    eigenincFound=no;
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    oldCPPFLAGS="$CPPFLAGS"
    CPPFLAGS="-I$EIGEN3_DIR $CPPFLAGS"
    AC_CHECK_HEADERS($EIGEN3_DIR/Eigen/Eigen, eigenincFound=yes)
    CPPFLAGS=$oldCPPFLAGS
    AC_LANG_RESTORE

    if (test $eigenincFound = no); then
      AC_MSG_RESULT(Eigen header files not found!)
      enableeigen=no;
    else
    	enableeigen=yes
    fi

    dnl If the Eigen headers were found, continue.
    if (test $enableeigen = yes); then
      EIGEN3_INCLUDES="-I$EIGEN3_DIR"
      CPPFLAGS="$EIGEN3_INCLUDES $CPPFLAGS"
      DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-eigen3=\"${EIGEN3_DIR}\""
      AC_DEFINE(HAVE_EIGEN, 1, [Flag indicating whether the library will be compiled with Eigen support])

      # Let us explicitly disable -Wshadow warnings since Eigen3 does not respect shadow declarations
      # and it results in an exorbitant amount of warnings during the build
      case "$CXXFLAGS" in
        *"shadow"*) CXXFLAGS="$CXXFLAGS -Wno-shadow"
      esac
    fi
  fi
  
  dnl Substitute the substitution variables
  AC_SUBST(EIGEN3_DIR) 
  AC_SUBST(EIGEN3_INCLUDES) 
  AC_SUBST(enableeigen)
])
