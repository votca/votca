HGVERSION:= $(shell hg -R $(top_builddir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_builddir)/.hg_archival.txt 2> /dev/null || true )

# Scripts generated at compile time
bin_SCRIPTS   = csg_inverse csg_call multi_g_rdf

if NO_RC_FILES
rcdir =
rcfiles =
else
rcdir = ${datadir}/votca/rc
rcfiles = csg.rc.csh csg.rc.bash
rc_DATA = $(rcfiles) csg-completion.bash
endif

# Extra files we need
EXTRA_DIST = csg.rc.csh.in csg.rc.bash.in csg_inverse.in csg_call.in multi_g_rdf.in csg-completion.bash

#hg version is determined at compiletime
csg_call: $(srcdir)/csg_call.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

csg_inverse: $(srcdir)/csg_inverse.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

multi_g_rdf: $(srcdir)/multi_g_rdf.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' $< > $@

csg.rc.csh: $(srcdir)/csg.rc.csh.in Makefile
	sed -e 's,@votcashare\@,${datadir}/votca,' $< > $@

csg.rc.bash: $(srcdir)/csg.rc.bash.in Makefile
	sed -e 's,@votcashare\@,${datadir}/votca,' $< > $@

CLEANFILES = $(bin_SCRIPTS) $(rcfiles) hgstamp *~

install-exec-hook:
	cd $(DESTDIR)$(bindir) && $(LN_S) -f multi_g_rdf multi_g_density

uninstall-hook:
	cd $(DESTDIR)$(bindir) && rm -f multi_g_density

hgstamp: dummy
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@

dummy: ;
