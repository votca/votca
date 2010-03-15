#! /bin/sh -e
aclocal -I config
autoheader
automake --add-missing --copy
autoconf 
