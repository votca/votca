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

#xml, gmx is needed because modules/io headers
libvotca_csg_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(GMX_CFLAGS) $(XML_CFLAGS) $(BOOST_CFLAGS)
#XML_LIBS come with VOTCA_TOOLS_LIBS
#boost is need in cgengine.h
libvotca_csg_la_LIBADD = modules/io/libio.la
libvotca_csg_la_LDFLAGS = -no-undefined $(VOTCA_TOOLS_LIBS) $(GMX_LIBS) $(BOOST_LIBS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_csg.pc

CLEANFILES = *~ hgstamp

#install-exec-hook:
#	rm -f $(DESTDIR)$(libdir)/libvotca_csg.la

version.o: hgstamp
version.lo: hgstamp
hgstamp: dummy 
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@
dummy: ;
