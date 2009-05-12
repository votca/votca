# die function
die () {
  echo "$*" 
  exit 1
}

# creating Makefile.am for share/scripts/inverse, this is automated here because we always forgot to update
dir=$PWD
cd share/scripts/inverse || die "error: cd shared/scripts/inverse failed"
rm -f Makefile.am
sed -ne '1,/%scripts/p' Makefile.am.in | sed -e '$d' > Makefile.am \
  || die "error creating share/scripts/inverse/Makefile.am"
ls *.sh *.pl csg_table *.m *.octave | sed -e 's/$/ \\/' -e '$s/ \\//' >> Makefile.am \
  || die "error creating share/scripts/inverse/Makefile.am"
sed -ne '/%scripts/,$p' Makefile.am.in | sed -e '1d' >> Makefile.am \
  || die "error creating share/scripts/inverse/Makefile.am"
cd $dir

# now do usual stuff
aclocal -I config
autoheader
automake --add-missing --copy
autoconf 

