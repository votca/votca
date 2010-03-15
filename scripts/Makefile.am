# Scripts generated at compile time
bin_SCRIPTS   = VOTCARC.csh VOTCARC.bash

# Extra files we need
EXTRA_DIST = VOTCARC.csh.in VOTCARC.bash.in

compdir = ${datadir}/votca
comp_DATA = completion.bash

VOTCARC.csh: $(srcdir)/VOTCARC.csh.in Makefile
	sed -e 's,@bindir\@,$(bindir),' \
	    -e 's,@libdir\@,$(libdir),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

VOTCARC.bash: $(srcdir)/VOTCARC.bash.in Makefile
	sed -e 's,@bindir\@,$(bindir),' \
	    -e 's,@libdir\@,$(libdir),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

CLEANFILES = $(bin_SCRIPTS) *~

