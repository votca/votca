# Scripts generated at compile time
bin_SCRIPTS   = CSGRC.csh CSGRC.bash csg_inverse csg_call

# Extra files we need
EXTRA_DIST = CSGRC.csh.in CSGRC.bash.in 

CSGRC.csh:  CSGRC.csh.in Makefile
	cat $(srcdir)/CSGRC.csh.in | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@votcashare\@,${datadir}/votca,' > CSGRC.csh

CSGRC.bash:  CSGRC.bash.in Makefile
	cat $(srcdir)/CSGRC.bash.in | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@votcashare\@,${datadir}/votca,' > CSGRC.bash

CLEANFILES = CSGRC.csh CSGRC.bash
