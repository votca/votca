HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_srcdir)/.hg_archival.txt 2> /dev/null || true )
AM_CPPFLAGS = -I$(top_srcdir)/include -DHGVERSION="\"$(HGVERSION)\""

SUBDIRS = calculators 

lib_LTLIBRARIES = libvotca_md2qm.la

libvotca_md2qm_la_SOURCES = \
	qmnblist.cc qmapplication.cc qmpair.cc  qmtopology.cc  statesaver.cc  version.cc calculatorfactory.cc

#gmx is needed because of headers from libio
libvotca_md2qm_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(VOTCA_CSG_CFLAGS) $(VOTCA_MOO_CFLAGS) $(VOTCA_KMC_CFLAGS)

libvotca_md2qm_la_LDFLAGS = -no-undefined $(VOTCA_TOOLS_LIBS) $(VOTCA_CSG_LIBS) $(VOTCA_MOO_LIBS) $(VOTCA_KMC_LIBS)
libvotca_md2qm_la_LIBADD = calculators/libcalculators.la

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_md2qm.pc

CLEANFILES = *~ hgstamp

install-exec-hook:
if NO_LA_FILES
	rm -f $(DESTDIR)$(libdir)/libvotca_md2qm.la
endif

version.o: hgstamp
version.lo: hgstamp
hgstamp: dummy 
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@
dummy: ;
