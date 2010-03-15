#! /bin/sh -e
aclocal -I config
autoheader
automake --add-missing --copy
autoconf 
[ -x ./src/libboost/bootstrap.sh ] && cd ./src/libboost && ./bootstrap.sh
