# Scripts generated at compile time
if NO_RC_FILES
bin_SCRIPTS =
else
bin_SCRIPTS   = VOTCARC.csh VOTCARC.bash
endif

# Extra files we need
EXTRA_DIST = VOTCARC.csh.in VOTCARC.bash.in

VOTCARC.csh: $(srcdir)/VOTCARC.csh.in Makefile
	sed -e 's,@bindir\@,$(bindir),' \
	    -e 's,@libdir\@,$(libdir),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

VOTCARC.bash: $(srcdir)/VOTCARC.bash.in Makefile
	sed -e 's,@bindir\@,$(bindir),' \
	    -e 's,@libdir\@,$(libdir),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

CLEANFILES = $(bin_SCRIPTS) *~

