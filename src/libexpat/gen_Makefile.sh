#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

find_to_make() {
  sort | \
  sed -e 's#^\./##' | \
  fmt | \
  sed -e 's#$# \\#' -e '$s# \\$##'
}

allow=" autom4te.cache config libs include "
for i in *; do
  [ -z "${allow//* $i *}" ] && continue 
  [ -d "$i" ] && die "Unknown dir $i"
done

echo Doing libs
cd libs
files=$(find . -type f -not -name "*.c" -and -not -name "Makefile*" -and -not -name "*.o" -and -not -name "*.l[oa]")
files=$(echo "$files" | grep -Eve '/.(libs|deps)/')

[ -n "$files" ] && die "Unknown files:\n$files"

rm -f Makefile.am
cat << eof >> Makefile.am
libvotca_expat_la_CPPFLAGS = -I\$(srcdir)/../include -DHAVE_EXPAT_CONFIG_H

libvotca_expat_la_LDFLAGS = -no-undefined

lib_LTLIBRARIES = libvotca_expat.la

libvotca_expat_la_SOURCES = \\
eof
files=$(find . -type f -name "*.c")
[ -z "${files}" ] && die "No cpp files found"
find . -type f -name "*.c" | find_to_make >> Makefile.am
echo >> Makefile.am
cat << eof >> Makefile.am
install-exec-hook:
if NO_LA_FILES
	rm -f \$(DESTDIR)\$(libdir)/libvotca_export.la
endif

eof

cd ..
echo Done with libs
echo
echo Doing headers
cd include
files=$(find . -type f -not -name "*.h" -and -not -name "Makefile*")
[ -n "$files" ] && die "Unknown files:\n$files"

rm -f Makefile.am
files=$(find . -type f -name "*.h")
[ -z "${files}" ] && die "No h files found"
echo -e "expatdir = \$(includedir)/votca/expat\n" >> Makefile.am
echo -e "expat_HEADERS = \\" >> Makefile.am
echo "$files" | find_to_make >> Makefile.am
cd ..
echo Done with headers
