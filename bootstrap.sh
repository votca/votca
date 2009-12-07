#!/bin/bash -e

# creating Makefile.am for share/scripts/inverse, this is automated here because we always forgot to update
make -C share/scripts/inverse -f Makefile.am.in Makefile.am

# now do usual stuff
aclocal -I config
autoheader
automake --add-missing --copy
autoconf 

