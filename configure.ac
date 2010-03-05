AC_INIT(votca-tools,1.0_rc4,bugs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])


AC_PROG_CXX
LT_INIT

AC_LANG_CPLUSPLUS

AC_HEADER_STDC
PKG_PROG_PKG_CONFIG

AC_ARG_WITH(votca_boost,
            [AC_HELP_STRING([--with-votca_boost],
                            [Build libvotca_boost, a build-in replacement for boost])],,
                             [with_votca_boost="no"])
# check for boost
if test "$with_votca_boost" = no; then
  AX_BOOST_BASE([1.33.1],[have_boost=yes],[have_boost=no])
fi

if test "$have_boost" = "no" -o "$want_boost" = "no"; then
  AC_MSG_NOTICE([
		
BOOST was not found/disabled: Building libvotca_boost
])
  with_votca_boost="yes"
fi

# check for fftw3
PKG_CHECK_MODULES([FFTW],fftw3)

# check for libxml2
PKG_CHECK_MODULES([XML],libxml-2.0)

# check for gsl
PKG_CHECK_MODULES([GSL],gsl)

AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
AM_CONDITIONAL(HAVE_TXT2TAGS,[test "$TXT2TAGS" != no])


if test "$with_votca_boost" = "yes"; then
  BOOST_CPPFLAGS='-I$$(top_srcdir)/src/libboost'
  AC_SUBST(BOOST_CPPFLAGS)
  BOOST_LDFLAGS="$$(top_srcdir)/src/libboost/libs/libvotca_boost.la"
  AC_SUBST(BOOST_LDFLAGS)
fi
AM_CONDITIONAL(VOTCA_BOOST, test "$with_votca_boost" = "yes")

AC_CONFIG_HEADERS(config.h)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libtools/Makefile])
AC_CONFIG_FILES([src/libtools/libvotca_tools.pc])

if test "$with_votca_boost" = "yes"; then
  AC_CONFIG_FILES([src/libboost/boost/algorithm/string/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/algorithm/string/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/algorithm/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/bind/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/config/abi/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/config/compiler/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/config/no_tr1/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/config/platform/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/config/stdlib/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/config/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/exception/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/exception/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/function/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/function/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/iterator/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/iterator/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/config/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/bcc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/bcc551/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/bcc_pre590/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/dmc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/gcc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/msvc60/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/msvc70/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/mwcw/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/no_ctps/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/no_ttp/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/plain/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessed/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/preprocessor/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/aux_/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/limits/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/vector/aux_/preprocessed/no_ctps/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/vector/aux_/preprocessed/plain/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/vector/aux_/preprocessed/typeof_based/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/vector/aux_/preprocessed/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/vector/aux_/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/vector/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/mpl/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/numeric/ublas/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/numeric/ublas/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/numeric/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/arithmetic/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/arithmetic/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/array/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/comparison/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/config/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/control/detail/dmc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/control/detail/edg/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/control/detail/msvc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/control/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/control/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/debug/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/detail/dmc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/facilities/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/iteration/detail/bounds/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/iteration/detail/iter/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/iteration/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/iteration/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/list/detail/dmc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/list/detail/edg/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/list/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/list/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/logical/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/punctuation/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/repetition/detail/dmc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/repetition/detail/edg/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/repetition/detail/msvc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/repetition/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/repetition/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/seq/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/seq/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/slot/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/slot/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/tuple/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/preprocessor/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/program_options/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/program_options/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/range/detail/vc6/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/range/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/range/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/serialization/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/serialization/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/smart_ptr/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/smart_ptr/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/type_traits/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/type_traits/msvc/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/type_traits/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/utility/Makefile])
  AC_CONFIG_FILES([src/libboost/boost/Makefile])
  AC_CONFIG_FILES([src/libboost/libs/detail/Makefile])
  AC_CONFIG_FILES([src/libboost/libs/program_options/Makefile])
  AC_CONFIG_FILES([src/libboost/libs/program_options/src/Makefile])
  AC_CONFIG_FILES([src/libboost/libs/smart_ptr/src/Makefile])
  AC_CONFIG_FILES([src/libboost/libs/smart_ptr/Makefile])
  AC_CONFIG_FILES([src/libboost/libs/Makefile])
  AC_CONFIG_FILES([src/libboost/Makefile])
  AC_CONFIG_FILES([src/libboost/libvotca_boost.pc])
fi

AC_OUTPUT
