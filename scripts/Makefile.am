HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_builddir)/.hg_archival.txt 2> /dev/null || true )

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
EXTRA_DIST = csg.rc.csh.in csg.rc.bash.in csg_inverse.in csg_call.in multi_g_rdf.in csg-completion.bash help2t2t.in votca-csg.t2t

#hg version is determined at compiletime
csg_call: $(srcdir)/csg_call.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@
	chmod +x $@

csg_inverse: $(srcdir)/csg_inverse.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' \
	    -e 's,@votcashare\@,${datadir}/votca,' $< > $@
	chmod +x $@

multi_g_rdf: $(srcdir)/multi_g_rdf.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' $< > $@
	chmod +x $@

multi_g_density: multi_g_rdf
	$(LN_S) -f $< $@

csg.rc.csh: $(srcdir)/csg.rc.csh.in Makefile
	sed -e 's,@votcashare\@,${datadir}/votca,' $< > $@

csg.rc.bash: $(srcdir)/csg.rc.bash.in Makefile
	sed -e 's,@votcashare\@,${datadir}/votca,' $< > $@

if HAVE_TXT2TAGS
man1_MANS = csg_call.man csg_inverse.man multi_g_rdf.man multi_g_density.man
man7_MANS = votca-csg.man

help2t2t: $(srcdir)/help2t2t.in Makefile hgstamp
	sed -e 's,@version\@,@VERSION@ $(HGVERSION),' $< > $@
	chmod +x $@

csg_inverse.t2t: csg_inverse help2t2t
	CSGSHARE=${top_builddir}/share ./help2t2t $< > $@

%.man: %.t2t
	${TXT2TAGS} -t man -i $< -o $@

%.t2t: % help2t2t
	./help2t2t $* > $@

else
man1_MANS=
man7_MANS=
endif

CLEANFILES = $(bin_SCRIPTS) $(rcfiles) $(man1_MANS) $(man7_MANS) hgstamp multi_g_density help2t2t csg_inverse.t2t *~

install-exec-hook:
	cd $(DESTDIR)$(bindir) && $(LN_S) -f multi_g_rdf multi_g_density

uninstall-hook:
	cd $(DESTDIR)$(bindir) && rm -f multi_g_density

BUILT_SOURCES = hgstamp

hgstamp: dummy
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@

dummy: ;
