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
files=$(find . -type f -not -name "*.cpp" -and -not -name "Makefile*" -and -not -name "*.Plo")
[ -n "$files" ] && die "Unknown files:\n$files"

rm -f Makefile.am
echo -e "AM_CPPFLAGS = -I\$(srcdir)/..\n" >>  Makefile.am
echo -e "lib_LTLIBRARIES = libvotca_boost.la\n" >> Makefile.am
echo -e "libvotca_boost_la_SOURCES = \\" >> Makefile.am
find . -type f -name "*.cpp" | grep -v detail | find_to_make >> Makefile.am
echo >> Makefile.am
echo -e "EXTRA_DIST = \\" >> Makefile.am
find . -type f -name "*.cpp" | grep detail | find_to_make >> Makefile.am
echo >> Makefile.am

cd ..
echo Done with libs
echo
echo Doing headers
cd boost
files=$(find . -type f -not -name "*.hpp" -a -not -name "Makefile*")
[ -n "$files" ] && die "Unknown files:\n$files"

rm -f Makefile.am
echo -e "boostincldir = \$(includedir)/votca/boost\n" >> Makefile.am
echo -e "nobase_boostincl_HEADERS = \\" >> Makefile.am
find . -type f -name "*.hpp" | find_to_make >> Makefile.am
echo >> Makefile.am

cd ..
echo Done with headers
