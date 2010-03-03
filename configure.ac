AC_INIT(votca-csg,1.0_rc4,bugs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_HEADERS(src/config.h)

AC_PROG_CXX
AC_LANG([C++])
LT_INIT
AC_HEADER_STDC
AC_PROG_LN_S
PKG_PROG_PKG_CONFIG

AX_GROMACS
AX_BOOST_BASE([1.33.1])
AX_BOOST_PROGRAM_OPTIONS([1.33.1])

AC_ARG_ENABLE(gmx4dev,
              [AC_HELP_STRING([--enable-gmx4dev],
              [Use gromacs 4 devel version])],, enable_gmx4dev=no)

if test "$enable_gmx4dev" = "yes"; then
  AC_DEFINE(GMX4DEV,,[Use gromacs 4 devel version])
fi

AC_ARG_ENABLE(all-static,
              [AC_HELP_STRING([--enable-all-static],
              [Make completely static binaries])],, enable_all_static=no)
# Dont add the -all-static flag until after the compiler test.

#look for votca tools
PKG_CHECK_MODULES([VOTCA_TOOLS],libvotca_tools)

#look for gsl
PKG_CHECK_MODULES([GSL],gsl)


if test "$enable_all_static" = "yes"; then
  AX_BROKEN_LIBXML
  LDFLAGS="$LDFLAGS -all-static"
fi

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
