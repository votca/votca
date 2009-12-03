HGVERSION:= $(shell grep -q node $(top_builddir)/hg_archival.txt 2> /dev/null || hg parents --template 'hgid: {node}' 2> /dev/null || true )

# Scripts generated at compile time
bin_SCRIPTS   = CSGRC.csh CSGRC.bash csg_inverse csg_call

# Extra files we need
EXTRA_DIST = CSGRC.csh.in CSGRC.bash.in csg_inverse.in csg_call.in

#hg version is determined at compiletime
csg_call: csg_call.in Makefile hgstamp
	cat $(srcdir)/csg_call.in | \
	sed 's,@version\@,@VERSION@ $(HGVERSION),' > csg_call

csg_inverse: csg_inverse.in Makefile hgstamp
	cat $(srcdir)/csg_inverse.in | \
	sed 's,@version\@,@VERSION@ $(HGVERSION),' > csg_inverse

CSGRC.csh:  CSGRC.csh.in Makefile hgstamp
	cat $(srcdir)/CSGRC.csh.in | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@votcashare\@,${datadir}/votca,' > CSGRC.csh

CSGRC.bash:  CSGRC.bash.in Makefile hgstamp
	cat $(srcdir)/CSGRC.bash.in | \
	sed 's,@bindir\@,$(bindir),' | \
	sed 's,@votcashare\@,${datadir}/votca,' > CSGRC.bash

CLEANFILES = CSGRC.csh CSGRC.bash csg_call csg_inverse

hgstamp: update_hgstamp
update_hgstamp:
	[ -f hgstamp ] || touch hgstamp
	echo $(HGVERSION) | cmp -s hgstamp - || echo $(HGVERSION) > hgstamp

