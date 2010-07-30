AC_DEFUN([AX_GROMACS], [
  AC_ARG_VAR([GMXLDLIB],[path to gromacs lib dir, usually set by "source GMXRC"])
  AC_ARG_WITH(libgmx,
    [AS_HELP_STRING([--with-libgmx@<:@=ARG@:>@],
      [use Gromacs(gmx) library, default single precision (ARG=libgmx/yes),
        but double precision (ARG=libgmx_d), mpi version (ARG=libgmx_mpi)
        or none (ARG=no or --without-libgmx) are possible as well @<:@ARG=yes@:>@])],,
  [with_libgmx=yes])
  if test "$with_libgmx" = "no"; then
    AC_SUBST(PKGGMX,"")
    AC_SUBST(PKGCFLAGSGMX,"")
    AC_SUBST(PKGLIBSGMX,"")
  else 
    if test -z "$withval" -o "$withval" = "yes"; then
      libgmx="gmx"
    else
      libgmx="${withval#lib}"
    fi
    PKG_CHECK_MODULES([GMX],[lib$libgmx],[:],[
      AC_MSG_CHECKING([GMXLDLIB])
      if test -z "$GMXLDLIB"; then
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([

Could not find GMXLDLIB environment variable, please source <gomacs-path>/bin/GMXRC 
or specify GMX_LIBS and GMX_CFLAGS
        ])
      else
        AC_MSG_RESULT([yes])
      fi
      AC_MSG_NOTICE([creating GMX_LIBS and GMX_CFLAGS from GMXLDLIB])
      if test -z "$GMX_LIBS"; then
        GMX_LIBS="-L$GMXLDLIB -l$libgmx"
        AC_MSG_NOTICE([setting GMX_LIBS   to "$GMX_LIBS"])
      else
        AC_MSG_NOTICE([GMX_LIBS was already set elsewhere to "$GMX_LIBS"])
      fi
    ])
    save_LIBS="$LIBS"
    save_CXX="$CXX"

    LIBS="$GMX_LIBS $LIBS"
    CXX="${SHELL-/bin/sh} ./libtool --mode=link $CXX"
    AC_MSG_CHECKING([for GromacsVersion in $GMX_LIBS])
    AC_TRY_LINK_FUNC(GromacsVersion,[AC_MSG_RESULT([yes])],[
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([

Could not link against GROMACS, please choose libraries to link against with --with-libgmx=XXX
or specify it in GMX_LIBS (e.g. export GMX_LIBS="-L<gromacs-path>/lib -lgmx")
or disable GROMACS support with --without-libgmx.

If you are using a mpi version of gromacs, make sure that CXX is something like mpic++.
(e.g. export CXX="mpic++" and export GMX_LIBS="-L<gromacs-path>/lib -lgmx_mpi")
      ])
    ])
    AC_MSG_CHECKING([for output_env_init in $GMX_LIBS])
    AC_TRY_LINK_FUNC(output_env_init,[
      AC_MSG_RESULT([yes])
      AC_MSG_NOTICE([

We are using a development version of gromacs, we hope you know what you are doing....
unexpected results ahead
      ])
      AC_DEFINE(GMX4DEV,1,[Use gromacs 4 devel version])
      gmxsub=""
      gmxheader="gromacs/tpxio.h"
    ],[
      AC_MSG_RESULT([no])
      gmxsub="/gromacs"
      gmxheader="tpxio.h"
    ])
    PKG_CHECK_EXISTS([lib$libgmx],[:],[
      if test -z "$GMX_CFLAGS"; then
        if test -z "${libgmx##*_d}"; then
          GMX_CFLAGS="-DGMX_DOUBLE -I$GMXLDLIB/../include$gmxsub"
	else
          GMX_CFLAGS="-I$GMXLDLIB/../include$gmxsub"
	fi
        AC_MSG_NOTICE([setting GMX_CFLAGS to "$GMX_CFLAGS"])
      else
        AC_MSG_NOTICE([GMX_CFLAGS was already set elsewhere to "$GMX_CFLAGS"])
      fi
    ])
    CXX="$save_CXX"
    LIBS="$save_LIBS"
    
    save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$GMX_CFLAGS $CPPFLAGS"
    AC_CHECK_HEADERS([$gmxheader],,[
      AC_MSG_ERROR([

Gromacs headers not found,
please make sure GMX_CFLAGS is pointing to <gomacs-path>/include for gromacs version >= 4.5
                                     or to <gomacs-path>/include/gromacs for gromacs version <= 4.0
      ])
    ])
  
    CPPFLAGS="$save_CPPFLAGS"
 
    dnl we need to do PKG_CHECK_EXISTS to know if libgmx pkg-config file
    dnl really exist, so that we can add it to our pkg-config files
    PKG_CHECK_EXISTS([lib$libgmx],[
      AC_SUBST(PKGGMX,"libgmx")
      AC_SUBST(PKGCFLAGSGMX,"")
      AC_SUBST(PKGLIBSGMX,"")
    ],[
      AC_SUBST(PKGGMX,"")
      AC_SUBST(PKGCFLAGSGMX,"$GMX_CFLAGS")
      AC_SUBST(PKGLIBSGMX,"$GMX_LIBS")
    ])
    AC_DEFINE(GMX,1,[Do we have gromacs])
  fi
  AM_CONDITIONAL(GMX,[test "$with_libgmx" != "no"])
])
