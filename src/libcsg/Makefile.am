HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_srcdir)/.hg_archival.txt 2> /dev/null || true )
AM_CPPFLAGS = -I$(top_srcdir)/include -DHGVERSION="\"$(HGVERSION)\""

SUBDIRS = modules 

lib_LTLIBRARIES = libvotca_csg.la

libvotca_csg_la_SOURCES = \
    beadlist.cc       exclusionlist.cc    nblist.cc         topology.cc          trajectorywriter.cc \
    cgengine.cc       neighbourlist.cc  topologymap.cc \
    cgmoleculedef.cc  map.cc              nematicorder.cc   topologyreader.cc \
    connectivity.cc   molecule.cc         numberdist.cc     trajectoryreader.cc \
    imcio.cc nblistgrid.cc version.cc 

#gmx is needed because of headers from libio
libvotca_csg_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(GMX_CFLAGS) $(BOOST_CFLAGS)
libvotca_csg_la_LIBADD = modules/io/libio.la

#boost_libs is needed in cgengine.cc from cgengine.h
#expat comes with VOTCA_TOOLS_LIBS
libvotca_csg_la_LDFLAGS = -no-undefined $(VOTCA_TOOLS_LIBS) $(GMX_LIBS) $(BOOST_LIBS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_csg.pc

CLEANFILES = *~ hgstamp

install-exec-hook:
if NO_LA_FILES
	rm -f $(DESTDIR)$(libdir)/libvotca_csg.la
endif

version.o: hgstamp
version.lo: hgstamp
hgstamp: dummy 
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@
dummy: ;
