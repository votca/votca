AM_CPPFLAGS = -I$(top_srcdir)/include

SUBDIRS = modules 

BUILT_SOURCES = hgversion.h

lib_LTLIBRARIES = libvotca_csg.la

libvotca_csg_la_SOURCES = hgversion.h \
beadlist.cc boundarycondition.cc cgengine.cc cgmoleculedef.cc \
csgapplication.cc exclusionlist.cc imcio.cc map.cc molecule.cc nblist.cc \
nblistgrid.cc nematicorder.cc openbox.cc orthorhombicbox.cc topology.cc \
topologymap.cc topologyreader.cc trajectoryreader.cc trajectorywriter.cc \
triclinicbox.cc version.cc

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
