AC_DEFUN([AX_BROKEN_LIBXML], [
  PKG_CHECK_MODULES([XML],libxml-2.0)
  AC_ARG_VAR([XML_STATIC_LIBS], [linker flags for static XML check, overriding pkg-config])
  if test -z "$XML_STATIC_LIBS"; then
    PKG_CHECK_EXISTS([libxml-2.0],
      [XML_STATIC_LIBS=`$PKG_CONFIG --libs --static libxml-2.0 2>/dev/null`],
      AC_MSG_FAILURE([Could not get libs for static linking of libxml])
    )
  fi
  AC_MSG_CHECKING([for broken libxml])
  CPPFLAGS="-static $XML_CFLAGS"
  dnl Sometimes libxml has a missing -lpthreads
  LIBS="$XML_LIBS $XML_STATIC_LIBS -lpthread"
  AC_RUN_IFELSE(
    AC_LANG_SOURCE([
#include <libxml/parser.h>

int main(int argc, char **argv) {
  xmlInitParser();
}
    ]),
    [AC_MSG_RESULT([no])],
    [
      AC_MSG_RESULT([yes])
      AC_MSG_FAILURE([Your libxml2 is broken for static linking, rebuild it without threads])]
  ,)
  CPPFLAGS="$save_CPPFLAGS"
  LIBS="$save_LIBS"
])
