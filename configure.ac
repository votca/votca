AC_INIT(votca-csg,1.0_rc5,bugs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.8 foreign])

AC_CONFIG_HEADERS(src/config.h)

AC_PROG_CXX
AC_LANG([C++])
LT_INIT
AC_HEADER_STDC
AC_PROG_LN_S
PKG_PROG_PKG_CONFIG

#libtool is needed by AX_GROMACS and AX_VOTCA_TOOLS
LT_OUTPUT

AC_ARG_ENABLE(all-static,
              [AC_HELP_STRING([--enable-all-static],
              [Make completely static binaries])],, enable_all_static=no)

if test "$enable_all_static" = "yes"; then
  BIN_LDFLAGS="-all-static $BIN_LDFLAGS"
  AC_SUBST(BIN_LDFLAGS)
fi

AC_ARG_ENABLE(la-files,
	[AC_HELP_STRING([--disable-la-files],[Do NOT install GNU libtool's la files])],,
	[enable_la_files=yes])

AM_CONDITIONAL(NO_LA_FILES,[test "$enable_la_files" = no])

AX_EXPAT

AX_VOTCA_BOOST

#look for votca tools
AX_VOTCA_TOOLS

AX_GROMACS

AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
AM_CONDITIONAL(HAVE_TXT2TAGS, test "$TXT2TAGS" != no)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/csg_boltzmann/Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/tools/Makefile])
AC_CONFIG_FILES([src/libcsg/Makefile])
AC_CONFIG_FILES([src/libcsg/libvotca_csg.pc])
AC_CONFIG_FILES([src/libcsg/modules/Makefile])
AC_CONFIG_FILES([src/libcsg/modules/io/Makefile])

AC_CONFIG_FILES([share/scripts/inverse/Makefile])
AC_CONFIG_FILES([share/scripts/Makefile])
AC_CONFIG_FILES([share/Makefile])
AC_CONFIG_FILES([scripts/Makefile])

AC_CONFIG_FILES([share/xml/Makefile])
AC_CONFIG_FILES([share/template/Makefile])

AC_OUTPUT
