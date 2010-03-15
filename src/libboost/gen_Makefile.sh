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

allow=" autom4te.cache libs boost config "
for i in *; do
  [ -z "${allow//* $i *}" ] && continue 
  [ -d "$i" ] && die "Unknown dir $i"
done

echo Doing libs
cd libs
files=$(find . -type f -not -name "*.cpp" -and -not -name "Makefile*" -and -not -name "*.o" -and -not -name "*.l[oa]")
files=$(echo "$files" | grep -Eve '/.(libs|deps)/')

[ -n "$files" ] && die "Unknown files:\n$files"

rm -f Makefile.am
cat << eof >> Makefile.am
libvotca_boost_la_CPPFLAGS = -I\$(srcdir)/..

libvotca_boost_la_LDFLAGS = -no-undefined

lib_LTLIBRARIES = libvotca_boost.la

libvotca_boost_la_SOURCES = \\
eof
files=$(find . -type f -name "*.cpp")
[ -z "${files}" ] && die "No cpp files found"
find . -type f -name "*.cpp" | grep -v detail | find_to_make >> Makefile.am
echo >> Makefile.am
echo -e "EXTRA_DIST = \\" >> Makefile.am
find . -type f -name "*.cpp" | grep detail | find_to_make >> Makefile.am
cat << eof >> Makefile.am

install-exec-hook:
if NO_LA_FILES
	rm -f \$(DESTDIR)\$(libdir)/libvotca_tools.la
endif

eof

cd ..
echo Done with libs
echo
echo Doing headers
cd boost
files=$(find . -type f -not -name "*.hpp" -and -not -name "Makefile*")
[ -n "$files" ] && die "Unknown files:\n$files"

rm -f Makefile.am
files=$(find . -type f -name "*.hpp")
[ -z "${files}" ] && die "No hpp files found"
blocks=$(( $(echo "$files" | wc -l) / 100 ))
echo Makeing $blocks blocks
for ((c=0;c<=$blocks;c++)); do
  echo -e "boost${c}dir = \$(includedir)/votca/boost\n" >> Makefile.am
  echo -e "nobase_boost${c}_HEADERS = \\" >> Makefile.am
  a=$((c*100+1))
  b=$((a+99))
  echo "$files" | sed -n "${a},${b}p" | find_to_make >> Makefile.am
  echo >> Makefile.am
done
cd ..
echo Done with headers
