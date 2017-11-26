dnl-----------------------------------------------------------------------------------
dnl HYPRE.M4
dnl   A helper for installing the custom hypre required by diablo
dnl-----------------------------------------------------------------------------------


dnl-----------------------------------------------------------------------------------
dnl CONFIGURE HYPRE
dnl   Untar the package, then run the necessary installation steps
dnl-----------------------------------------------------------------------------------
AC_DEFUN([CONFIGURE_HYPRE],[
  # Call package Download/Configure/Installation procedures for HYPRE, if requested by user
  PPREFIX=HYPRE

  # CLI option for HYPRE
  AC_MSG_CHECKING([if HYPRE support is enabled])
  AC_ARG_WITH(hypre, 
  [AS_HELP_STRING([--with-hypre@<:@=DIR@:>@], [Specify HYPRE library to use for parallel solvers])
  AS_HELP_STRING([--without-hypre], [Disable support for native HYPRE solvers])],
  [if (test "x$withval" != "x" && test "x$withval" != "xno"); then
    HYPRE_DIR=$withval
    enablehypre=yes
    DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-hypre=\"${withval}\""
  fi
  ], [HYPRE_DIR=$HYPRE_DIR])
  if (test "x" != "x$HYPRE_DIR" && test "xno" != "x$HYPRE_DIR"); then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    # Reset the directory since we do not want to configure HYPRE
    HYPRE_DIR=""
    enablehypre=no
  fi

  # Supported Hypre versions: 2.13.0, 2.11.2, 2.10.1
  # Arguments: 1) Default Version Number, 2) Download by default ?
  AUSCM_CONFIGURE_DOWNLOAD_HYPRE([2.13.0],[no])

  # Set some variables
  if (test "x$enablehypre" != "xno"); then
    HYPRE_LIBS="-L$HYPRE_DIR/lib -lHYPRE"
    HYPRE_INCLUDES="-I$HYPRE_DIR/include"

    tmpLIBS=$LIBS
    AC_LANG_SAVE
    AC_LANG_PUSH([C])
    LIBS="$MOAB_LDFLAGS $HYPRE_LIBS $MOAB_LIBS"
    # Check if it worked
    PREFIX_PRINT([Verifying libHYPRE in -L$HYPRE_DIR/lib...])
    AC_CHECK_LIB([HYPRE],[main],
        [MOAB_LDFLAGS="$MOAB_LDFLAGS -L$HYPRE_DIR/lib"
        MOAB_LIBS="$HYPRE_LIBS $MOAB_LIBS"
        enablehypre=yes],
        [enablehypre=no],
        [$MOAB_LIBS]
      )
    LIBS=$tmpLIBS
    AC_LANG_POP([C])
  fi
 
  # Return some variables
  AC_SUBST(enablehypre)
  AC_SUBST(HYPRE_DIR)
  AC_SUBST(HYPRE_INCLUDES)
  AC_SUBST(HYPRE_LIBS)
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED SETUP PREPROCESS HYPRE
dnl   Figure out what needs to be done to get a valid HYPRE installation.
dnl   Arguments: [PACKAGE, SRC_DIR, INSTALL_DIR, NEED_CONFIGURATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_SETUP_PREPROCESS_HYPRE],[
  # uncompress and configure PACKAGE
  hypre_src_dir="$2/src"
  hypre_build_dir="$2/src"
  hypre_install_dir="$3"
  hypre_archive_name="$4"
  PPREFIX=HYPRE

  if (test ! -d "$hypre_src_dir" || test ! -f "$hypre_src_dir/configure" ); then
  	AC_MSG_ERROR([Invalid source configuration for HYPRE. Source directory $hypre_src_dir is invalid])
  fi

  # determine what steps we need to take to install hypre
  hypre_configured=false
  hypre_made=false
  hypre_installed=false
  if (test ! -d "$hypre_build_dir" ); then
    AS_MKDIR_P( $hypre_build_dir )
  else
    if (test -f "$hypre_build_dir/HYPRE_config.h" ); then
      hypre_configured=true
      if (test -f "$hypre_build_dir/lib/libHYPRE.a" || test -f "$hypre_build_dir/lib/libHYPRE.so" || test -f "$hypre_build_dir/lib/libHYPRE.dylib" ); then
        hypre_made=true
        if (test -f "$hypre_install_dir/include/HYPRE.h"); then
          hypre_installed=true
        fi
      fi
    fi
  fi
  # send the information back
  AS_IF([ ! $hypre_configured || $need_configuration ], [need_configuration=true], [need_configuration=false])
  AS_IF([ ! $hypre_made || $need_configuration ], [need_build=true], [need_build=false])
  AS_IF([ ! $hypre_installed || $need_configuration ], [need_installation=true], [need_installation=false])
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED SETUP POSTPROCESS HYPRE
dnl   Dummy macro to fit standard call pattern.  Tells MOAB we have HYPRE.
dnl   Arguments: [PACKAGE, SRC_DIR, INSTALL_DIR, NEED_CONFIGURATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_SETUP_POSTPROCESS_HYPRE],[
  # we have already checked configure/build/install logs for errors before getting here..
  enablehypre=yes
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED CONFIGURE HYPRE
dnl   Runs configure for HYPRE and looks for header files.
dnl   Arguments: [NEED_CONFIGURATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_CONFIGURE_HYPRE],[
if [ $1 ]; then
  # configure HYPRE
  if [ $need_configuration ]; then
    # configure PACKAGE with a minimal build: MPI
    export CC=$CC CXX=$CXX FC=$FC F77=$F77 CFLAGS="$CFLAGS -fPIC -DPIC" CXXFLAGS="$CXXFLAGS -fPIC -DPIC" FCFLAGS="$FCFLAGS -fPIC" FFLAGS="$FFLAGS -fPIC" LDFLAGS="$LDFLAGS"
    configure_command="$hypre_src_dir/configure --prefix=$hypre_install_dir --libdir=$hypre_install_dir/lib --with-pic=1 --enable-cxx"
    if (test "$enabledebug" != "no"); then
      configure_command="$configure_command --enable-debug"
    fi
    if (test "$enablefortran" != "no"); then
      configure_command="$configure_command --enable-fortran"
    fi
    if (test "$enable_shared" != "no"); then
      configure_command="$configure_command --enable-shared"
    fi
      
    if (test "x$MPI_DIR" != "x"); then
      configure_command="$configure_command --with-mpi=$MPI_DIR"
    fi  
    if (test "x$BLAS_DIR" != "x"); then
      configure_command="$configure_command --with-blas-lib-dirs=$BLAS_DIR --with-blas-libs=\"-lblas\""
    fi
    if (test "x$LAPACK_DIR" != "x"); then
      configure_command="$configure_command --with-lapack-lib-dirs=$LAPACK_DIR --with-lapack-libs=\"-llapack\""
    fi

    eval "echo 'Using configure command :==> cd $hypre_build_dir ; $configure_command > $hypre_src_dir/../../config_hypre.log ; cd -' > $hypre_src_dir/../../config_hypre.log"
    ##echo "Trying to use configure command:==> cd $hypre_build_dir ; $configure_command > $hypre_src_dir/../config_hypre.log ; cd -"
    PREFIX_PRINT(Configuring with default options  {debug=$enabledebug fortran=$enablefortran shared=$enable_shared BLAS=$BLAS_DIR LAPACK=$LAPACK_DIR} )
    eval "cd $hypre_build_dir ; $configure_command >> $hypre_src_dir/../config_hypre.log 2>&1 ; cd $MOAB_BUILD_DIR"
  fi

  # check if configuration - current or previous was successful
  if (test ! -f "$hypre_build_dir/HYPRE_config.h" ); then
	  AC_MSG_ERROR([HYPRE configuration was unsuccessful. Please refer to $hypre_build_dir/config.log and $hypre_src_dir/../../config_hypre.log for further details.])
  fi
  hypre_configured=true
fi
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED BUILD HYPRE
dnl   Calls make on HYPRE and looks for libraries.
dnl   Arguments: [NEED_BUILD)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_BUILD_HYPRE],
[
  # if we need to build then call make all
  if [ $1 ]; then
    if [ $recompile_and_install || $need_build ]; then
      PREFIX_PRINT(Building the sources in parallel )
      hypre_makelog="`cd $hypre_build_dir ; make all -j4 > $hypre_src_dir/../../make_hypre.log 2>&1 ; cd $MOAB_BUILD_DIR`"
    fi
  fi
  # check if it worked
  if (test -f "$hypre_build_dir/lib/libHYPRE.a" || test -f "$hypre_build_dir/lib/libHYPRE.so" || test -f "$hypre_build_dir/lib/libHYPRE.dylib") ; then
    hypre_made=true
  else
    AC_MSG_ERROR([HYPRE build was unsuccessful. Please refer to $hypre_src_dir/../../make_hypre.log for further details.])
  fi
])


dnl ---------------------------------------------------------------------------
dnl AUTOMATED INSTALL HYPRE
dnl   Calls make install on HYPRE and checks for libhypre.settings
dnl   Arguments: [NEED_INSTALLATION)
dnl ---------------------------------------------------------------------------
AC_DEFUN([AUTOMATED_INSTALL_HYPRE],
[
  # if we need to install then call make install
  if [ $1 ]; then
    if [ $recompile_and_install ]; then
      if [ $hypre_installed ]; then
        hypre_installlog="`cd $hypre_build_dir ; make uninstall > $hypre_src_dir/../../uninstall_hypre.log 2>&1 ; cd $MOAB_BUILD_DIR`"
      fi
      PREFIX_PRINT(Installing the headers and libraries in to directory ($hypre_install_dir) )
      hypre_installlog="`cd $hypre_build_dir ; make install > $hypre_src_dir/../../install_hypre.log 2>&1 ; cd $MOAB_BUILD_DIR`"
    fi
  fi
  # check if it worked
  if (test -f "$hypre_install_dir/include/HYPRE.h"); then
    hypre_installed=true
  else
	  AC_MSG_ERROR([HYPRE installation was unsuccessful. Please refer to $hypre_src_dir/../../install_hypre.log for further details.])
  fi
])

