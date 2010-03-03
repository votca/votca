HGVERSION:= $(shell hg -R $(top_builddir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_builddir)/.hg_archival.txt 2> /dev/null || true )

# Scripts generated at compile time
bin_SCRIPTS   = CSGRC.csh CSGRC.bash csg_inverse csg_call multi_g_rdf

comp_DATA = completion.bash
compdir = ${datadir}/votca

# Extra files we need
EXTRA_DIST = CSGRC.csh.in CSGRC.bash.in csg_inverse.in csg_call.in multi_g_rdf.in $(comp_DATA)

#hg version is determined at compiletime
csg_call: $(srcdir)/csg_call.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' $< > $@

csg_inverse: $(srcdir)/csg_inverse.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' $< > $@

multi_g_rdf: $(srcdir)/multi_g_rdf.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' $< > $@

CSGRC.csh: $(srcdir)/CSGRC.csh.in Makefile hgstamp
	sed -e 's,@bindir\@,$(bindir),' \
	    -e 's,@libdir\@,$(libdir),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

CSGRC.bash: $(srcdir)/CSGRC.bash.in Makefile hgstamp
	sed -e 's,@bindir\@,$(bindir),' \
	    -e 's,@libdir\@,$(libdir),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@

CLEANFILES = CSGRC.csh CSGRC.bash csg_call csg_inverse multi_g_rdf *~

install-exec-hook:
	$(LN_S) multi_g_rdf $(DESTDIR)$(bindir)/multi_g_density

hgstamp: dummy 
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@
dummy: ;
