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

    dnl Try to find pkg-files even if PKG_CONFIG_PATH was not set
    save_PKG_CONFIG_PATH="$PKG_CONFIG_PATH"
    AC_MSG_CHECKING([for GMX pkg-config files])
    PKG_CHECK_EXISTS([lib$libgmx],[
      AC_MSG_RESULT([yes])
    ],[
      AC_MSG_RESULT([no])
      AC_MSG_CHECKING([for GMXLDLIB])
      if test -z "$GMXLDLIB"; then
        AC_MSG_RESULT([no])
        AC_MSG_WARN([Could not find GMXLDLIB environment variable, assuming you have gromacs files in defaults paths])
      else
        AC_MSG_RESULT([yes])
      fi
      
      dnl in the case gromacs has pkg-config file, but PKG_CONFIG_PATH not set
      if test -n "$GMXLDLIB"; then
        AC_MSG_CHECKING([for GMX pkg-config dir])
        if test -d "$GMXLDLIB/pkgconfig"; then
          AC_MSG_RESULT([yes])
          export PKG_CONFIG_PATH="$GMXLDLIB/pkgconfig:$PKG_CONFIG_PATH"
        else
          AC_MSG_RESULT([no])
        fi
      fi
    ])

    dnl Again PKG_CHECK_MODULES after we added $GMXLDLIB/pkgconfig to PKG_CONFIG_PATH
    PKG_CHECK_MODULES([GMX],[lib$libgmx],[:],[
      dnl Do not overwrite user settings
      if test -z "$GMX_LIBS"; then
        if test -n "$GMXLDLIB"; then
          AC_MSG_NOTICE([creating GMX_LIBS from GMXLDLIB])
          GMX_LIBS="-L$GMXLDLIB -l$libgmx"
        else
          GMX_LIBS="-l$libgmx"
        fi
        AC_MSG_NOTICE([setting GMX_LIBS to "$GMX_LIBS"])
      fi
    ])

    save_LIBS="$LIBS"
    save_CXX="$CXX"
    LIBS="$GMX_LIBS $LIBS"
    CXX="${SHELL-/bin/sh} ./libtool --mode=link $CXX"
    AC_MSG_CHECKING([for GromacsVersion in $GMX_LIBS])
    AC_TRY_LINK_FUNC(GromacsVersion,[
      AC_MSG_RESULT([yes])
      GMX_VERSION="33"
    ],[
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([

Could not link against GROMACS, please choose libraries to link against with --with-libgmx=XXX
or specify it in GMX_LIBS (e.g. export GMX_LIBS="-L<gromacs-path>/lib -lgmx")
or disable GROMACS support with --without-libgmx.

If you are using a mpi version of gromacs, make sure that CXX is something like mpic++.
(e.g. export CXX="mpic++" and export GMX_LIBS="-L<gromacs-path>/lib -lgmx_mpi")
      ])
    ])
    AC_MSG_CHECKING([for init_mtop in $GMX_LIBS])
    AC_TRY_LINK_FUNC(init_mtop,[
      AC_MSG_RESULT([yes])
      GMX_VERSION="40"
    ],[
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([

Your version of GROMACS is too old, please update at least to version 4.0
      ])
    ])
    AC_MSG_CHECKING([for output_env_init in $GMX_LIBS])
    AC_TRY_LINK_FUNC(output_env_init,[
      AC_MSG_RESULT([yes])
      GMX_VERSION="45"
      gmxheader="gromacs/tpxio.h"
    ],[
      AC_MSG_RESULT([no])
      if test -n "$GMXLDLIB"; then
        gmxincldir="-I$GMXLDLIB/../include/gromacs"
      else
        gmxincldir="-I/usr/include/gromacs"
      fi
      gmxheader="tpxio.h"
    ])

dnl     enable this if gromacs 5.0 development starts
dnl     AC_MSG_CHECKING([for output_env_init in $GMX_LIBS])
dnl     AC_TRY_LINK_FUNC(output_env_init,[
dnl       AC_MSG_RESULT([yes])
dnl       AC_MSG_NOTICE([
dnl
dnl  We are using a development version of gromacs, we hope you know what you are doing....
dnl unexpected results ahead
dnl       ])
dnl       GMX_VERSION="50"
dnl     ],[
dnl       AC_MSG_RESULT([no])
dnl     ])
    CXX="$save_CXX"
    LIBS="$save_LIBS"

    PKG_CHECK_EXISTS([lib$libgmx],[:],[
      dnl Do not overwrite user settings
      if test -z "$GMX_CFLAGS"; then
        if test -n "$GMXLDLIB"; then
          AC_MSG_NOTICE([creating GMX_CFLAGS from GMXLDLIB])
        fi
        GMX_CFLAGS="$gmxincldir"
        dnl in case user wants double precision
        if test -z "${libgmx##*_d}"; then
          GMX_CFLAGS="-DGMX_DOUBLE $GMX_CFLAGS"
        fi
        AC_MSG_NOTICE([setting GMX_CFLAGS to "$GMX_CFLAGS"])
      fi
    ])
    
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
 
    dnl restore old PKG_CONFIG_PATH, otherwise pkg-config will search for libgmx.pc
    dnl when building against libvotca_csg and will not find it again
    PKG_CONFIG_PATH="$save_PKG_CONFIG_PATH"
    
    dnl we need to do PKG_CHECK_EXISTS to know if libgmx pkg-config file
    dnl really exist, so that we can add it to our pkg-config files
    PKG_CHECK_EXISTS([lib$libgmx],[
      AC_SUBST(PKGGMX,"lib$libgmx")
      AC_SUBST(PKGCFLAGSGMX,"")
      AC_SUBST(PKGLIBSGMX,"")
    ],[
      AC_SUBST(PKGGMX,"")
      AC_SUBST(PKGCFLAGSGMX,"$GMX_CFLAGS")
      AC_SUBST(PKGLIBSGMX,"$GMX_LIBS")
    ])
    AC_DEFINE_UNQUOTED(GMX,[$GMX_VERSION],[Used gromacs version])
  fi
  AM_CONDITIONAL(GMX,[test "$with_libgmx" != "no"])
])
