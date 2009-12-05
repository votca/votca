AC_INIT(votca-tools,1.0_rc1,devs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])


AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG_CPLUSPLUS

AC_HEADER_STDC

# check for boost
AX_BOOST_BASE([1.33.1])

# check for fftw3
AX_FFTW3

# check for libxml2
ACI_PACKAGE_SAVED(XML,[XML library (version XML2)],
                  xml2-config,,,
                  AC_MSG_ERROR([Cannot find the XML2 libs. Is xml2-config in your path ?])
                 )

AX_PATH_GSL(,,AC_MSG_ERROR([Cannot find GSL. Please make sure it is installed and check your flags]))

# check for gromacs
# AX_GROMACS_HEADERS

AC_PROG_MAKE_SET

AC_CONFIG_HEADERS(config.h)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libtools/Makefile])

AC_OUTPUT
