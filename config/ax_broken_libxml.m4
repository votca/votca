AC_DEFUN([AX_BROKEN_LIBXML], [
  if test -z "$XML_LIBS"; then
    PKG_CHECK_MODULES([XML],libxml-2.0)
    PKG_CHECK_EXISTS([libxml-2.0],
      [XML_LIBS=`$PKG_CONFIG --libs --static libxml-2.0 2>/dev/null`],
      AC_MSG_FAILURE([Could not get libs for static linking of libxml])
    )
  fi
  AC_MSG_CHECKING([for broken libxml])
  CPPFLAGS="-static $XML_CFLAGS"
  dnl Sometimes libxml has a missing -lpthreads
  LIBS="$2 $XML_LIBS -lpthread"
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
