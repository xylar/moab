# ===========================================================================
#          http://www.gnu.org/software/autoconf-archive/ax_blas.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the BLAS linear-algebra
#   interface (see http://www.netlib.org/blas/). On success, it sets the
#   BLAS_LIBS output variable to hold the requisite library linkages.
#
#   To link with BLAS, you should link with:
#
#     $BLAS_LIBS $LIBS $FLIBS
#
#   in that order. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_BLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   Many libraries are searched for, from ATLAS to CXML to ESSL. The user
#   may also use --with-blas=<lib> in order to use some specific BLAS
#   library <lib>. In order to link successfully, however, be aware that you
#   will probably need to use the same Fortran compiler (which can be set
#   via the F77 env. var.) as was used to compile the BLAS library.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a BLAS library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_BLAS.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 11

AU_ALIAS([ACX_BLAS], [AX_BLAS])
AC_DEFUN([AX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
ax_blas_ok=no

AC_ARG_WITH(blas,
  [AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
  yes | "") ;;
  no) ax_blas_ok=disable ;;
  -* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
  *) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_FC_FUNC(sgemm)
AC_FC_FUNC(dgemm)

ax_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $ax_blas_ok = no; then
if test "x$BLAS_LIBS" != "x"; then
  save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
  AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes], [BLAS_LIBS=""])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS"
  AC_MSG_CHECKING([if $sgemm is being linked in already])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-latlas $LIBS"
  AC_MSG_CHECKING([for $sgemm and $dgemm in -latlas])
  AC_TRY_LINK_FUNC(ATL_xerbla, 
          [LIBS="-lf77blas $LIBS";
           AC_TRY_LINK_FUNC($sgemm,
               [LIBS="-lcblas $LIBS";
                AC_TRY_LINK_FUNC(cblas_dgemm,
                   [ax_blas_ok=yes;
                    BLAS_LIBS="-lcblas -lf77blas -latlas"
                   ]) 
                ]) 
          ])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-lblas $LIBS"
  AC_MSG_CHECKING([for $sgemm and $dgemm in PhiPACK libraries])
  AC_TRY_LINK_FUNC($sgemm, 
          [LIBS="-ldgemm $LIBS";
           AC_TRY_LINK_FUNC($dgemm,
               [LIBS="-lsgemm $LIBS";
                AC_TRY_LINK_FUNC($sdgemm,
                   [ax_blas_ok=yes;
                    BLAS_LIBS="-lsgemm -ldgemm -lblas"
                   ]) 
                ]) 
          ])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in Intel MKL library?
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-lmkl $LIBS"
  AC_MSG_CHECKING([for $sgemm in -lmkl])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lmkl"])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in Apple vecLib library?
if (test "x$target_vendor" = "xapple" && test "x$BLAS_LIBS" = "x"); then

  if test $ax_blas_ok = no; then
    save_LIBS="$LIBS"; LIBS="-framework vecLib $LIBS"
    AC_MSG_CHECKING([for $sgemm in -framework vecLib])
    AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-framework vecLib"])
    AC_MSG_RESULT($ax_blas_ok)
    LIBS="$save_LIBS"
  fi
  if test $ax_blas_ok = no; then
    save_LIBS="$LIBS"; LIBS="-framework Accelerate $LIBS"
    AC_MSG_CHECKING([for $sgemm in -framework Accelerate])
    AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-framework Accelerate"])
    AC_MSG_RESULT($ax_blas_ok)
    LIBS="$save_LIBS"
  fi

fi

# BLAS in Alpha CXML library?
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-lcxml $LIBS"
  AC_MSG_CHECKING([for $sgemm in -lcxml])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lcxml"])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-ldxml $LIBS"
  AC_MSG_CHECKING([for $sgemm in -ldxml])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-ldxml"])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in Sun Performance library?
if test $ax_blas_ok = no; then
  if test "x$GCC" != xyes; then # only works with Sun CC
    AC_CHECK_LIB(sunmath, acosp,
      [AC_CHECK_LIB(sunperf, $sgemm,
        [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 ax_blas_ok=yes],[],[-lsunmath])])
  fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-lscs $LIBS"
  AC_MSG_CHECKING([for $sgemm in -lscs])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lscs"])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in SGIMATH library?
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-lcomplib.sgimath $LIBS"
  AC_MSG_CHECKING([for $sgemm in -lcomplib.sgimath])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lcomplib.sgimath"])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
  AC_CHECK_LIB(blas, $sgemm,
    [AC_CHECK_LIB(essl, $sgemm,
      [ax_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
      [], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $ax_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="-lblas $LIBS"
  AC_MSG_CHECKING([for $sgemm in -lblas])
  AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lblas"])
  AC_MSG_RESULT($ax_blas_ok)
  LIBS="$save_LIBS"
fi

AC_SUBST(BLAS_LIBS)

LIBS="$ax_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_blas_ok" = xyes; then
  AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.])
  AC_MSG_NOTICE([Found BLAS library])
  $1
else
  ax_blas_ok=no
  AC_MSG_ERROR([BLAS library not found])
  $2
fi
])dnl AX_BLAS
