AM_CPPFLAGS = -I$(top_srcdir)/include

SUBDIRS = modules 

BUILT_SOURCES = hgversion.h

lib_LTLIBRARIES = libvotca_csg.la

libvotca_csg_la_SOURCES = hgversion.h \
    beadlist.cc       exclusionlist.cc    nblist.cc         topology.cc          trajectorywriter.cc \
    cgengine.cc       topologymap.cc \
    cgmoleculedef.cc  map.cc              nematicorder.cc   topologyreader.cc \
    molecule.cc         trajectoryreader.cc \
    imcio.cc nblistgrid.cc version.cc csgapplication.cc

#gmx is needed because of headers from libio and in version.cc
libvotca_csg_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(GMX_CFLAGS)
libvotca_csg_la_LIBADD = modules/io/libio.la

libvotca_csg_la_LDFLAGS = -no-undefined $(VOTCA_TOOLS_LIBS) $(GMX_LIBS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_csg.pc

CLEANFILES = *~ hgversion.h

install-exec-hook:
if NO_LA_FILES
	rm -f $(DESTDIR)$(libdir)/libvotca_csg.la
endif

HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_srcdir)/.hg_archival.txt 2> /dev/null || true )
hgversion.h: dummy
	[ -f $@ ] || touch $@
	echo 'static const std::string hgversion = "$(HGVERSION)";' | cmp -s $@ - \
	  || echo 'static const std::string hgversion = "$(HGVERSION)";' > $@

.PHONY: dummy
dummy: ;
