AC_INIT(votca-csg,1.0_rc1,devs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_HEADERS(src/config.h)

AC_PROG_CXX
AC_LANG([C++])
AC_PROG_RANLIB
AC_HEADER_STDC
AC_PROG_MAKE_SET

# check for fftw3
AX_FFTW3

# check for libxml2
ACI_PACKAGE_SAVED(XML,[XML library (version XML2)],
                  xml2-config,,,
                  AC_MSG_WARN(Cannot find the XML2 libs. Is xml2-config in your path ?)
                 )

AX_PATH_GSL(,,AC_MSG_ERROR([Cannot find GSL. Please make sure it is installed and check your flags]))
AX_GROMACS
AX_BOOST([1.33.1])

AX_VOTCA_TOOLS
#check for boost
#AX_BOOST_PROGRAM_OPTIONS


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/csg_boltzmann/Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/tools/Makefile])
AC_CONFIG_FILES([src/libcsg/Makefile])
AC_CONFIG_FILES([src/libcsg/version.cc])

AC_CONFIG_FILES([share/scripts/inverse/Makefile])
AC_CONFIG_FILES([share/scripts/Makefile])
AC_CONFIG_FILES([share/Makefile])
AC_CONFIG_FILES([scripts/Makefile])

AC_CONFIG_FILES([share/xml/Makefile])
AC_CONFIG_FILES([share/template/Makefile])


AC_OUTPUT
