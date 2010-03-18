AC_DEFUN([AX_GROMACS], [
  AC_ARG_WITH(libgmx,[
    AS_HELP_STRING([--with-libgmx=NAME],
      [Compile against this Gromacs(gmx) library, default libgmx, but libgmx_d or libgmx_mpi is possible as well.
    ])
  ],,[with_libgmx=yes])
  if test "$with_libgmx" = "no"; then
    AC_MSG_ERROR([You cannot compile without libgmx])
  fi
  if test -z "$withval" -o "$withval" = "yes"; then
    libgmx="gmx"
  else
    libgmx="${withval#lib}"
  fi
  PKG_CHECK_MODULES([GMX],[lib$libgmx],[:],[
    AC_ARG_VAR([GMXLDLIB],[path to gromacs lib dir, usually set by "source GMXRC"])
    AC_MSG_CHECKING([GMXLDLIB])
    if test -z "$GMXLDLIB"; then
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([

Could not find GMXLDLIB environment variable, please source <gomacs-path>/bin/GMXRC 
or specify GMX_LIBS and GMX_CLFAGS
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
    if test -z "$GMX_CFLAGS"; then
      GMX_CFLAGS="-I$GMXLDLIB/../include/gromacs"
      AC_MSG_NOTICE([setting GMX_CFLAGS to "$GMX_CFLAGS"])
    else
      AC_MSG_NOTICE([GMX_CFLAGS was already set elsewhere to "$GMX_CFLAGS"])
    fi
  ])
  save_CPPFLAGS="$CPPFLAGS"
  save_LIBS="$LIBS"

  CPPFLAGS="$GMX_CFLAGS $CPPFLAGS"
  LIBS="$GMX_LIBS $LIBS"
  AC_CHECK_HEADERS([tpxio.h],,[
    AC_MSG_ERROR([

Gromacs headers not found,
please make sure GMX_CFLAGS is pointing to <gomacs-path>/include/gromacs 
Do not forget the /gromacs due to bug in gromacs headers!
    ])
  ])
  AC_MSG_CHECKING([for GromacsVersion in $GMX_LIBS])
  AC_TRY_LINK_FUNC(GromacsVersion,[AC_MSG_RESULT([yes])],[
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([

Could not link against GROMACS,please choose libraries to link against with --with-libgmx=XXX
or specify it in GMX_LIBS (e.g. export GMX_LIBS="-L<gromacs-path>/lib -lgmx").

If you are using a mpi version of gromacs, make sure that CXX is something like mpic++.
(e.g. export CXX="mpic++" and export GMX_LIBS="-L<gromacs-path>/lib -lgmx_mpi")
    ])
  ])
  AC_CHECK_HEADERS([types/oenv.h],[
    AC_MSG_CHECKING([for output_env_init in $GMX_LIBS])
    AC_TRY_LINK_FUNC(output_env_init,[AC_MSG_RESULT([yes])],[
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([

We found oenv.h header, which exist only in gromacs 4 devel, but we did find
output_env_init libgmx
      ])
    ])
    AC_MSG_NOTICE([

We are using a development version of gromacs, we hope you know what you are doing....
unexpected results ahead
    ])
    AC_DEFINE(GMX4DEV,,[Use gromacs 4 devel version])
  ])
  CPPFLAGS="$save_CPPFLAGS"
  LIBS="$save_LIBS"
  PKG_CHECK_EXISTS([lib$libgmx],[
    PKGGMX="libgmx"
    PKGCFLAGSGMX=""
    PKGLIBSGMX=""
  ],[
    PKGGMX=""
    PKGCFLAGSGMX="$GMX_CFLAGS"
    PKGLIBSGMX="$GMX_LIBS"
  ])
  AC_SUBST(PKGGMX)
  AC_SUBST(PKGCFLAGSGMX)
  AC_SUBST(PKGLIBSGMX)
])
