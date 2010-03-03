HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_srcdir)/.hg_archival.txt 2> /dev/null || true )
AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(GMX_CFLAGS) -I$(top_srcdir)/include -DHGVERSION="\"$(HGVERSION)\""

SUBDIRS = modules 

lib_LTLIBRARIES = libvotca_csg.la

libvotca_csg_la_SOURCES = \
    beadlist.cc       exclusionlist.cc    nblist.cc         topology.cc          trajectorywriter.cc \
   cgengine.cc       neighbourlist.cc  topologymap.cc \
   cgmoleculedef.cc  map.cc              nematicorder.cc   topologyreader.cc \
   connectivity.cc   molecule.cc         numberdist.cc     trajectoryreader.cc \
   imcio.cc nblistgrid.cc version.cc \
   modules/io/gmxtopologyreader.cc    modules/io/grotopologyreader.cc  modules/io/xmltopologyreader.cc \
    modules/io/gmxtrajectoryreader.cc  modules/io/growriter.cc          modules/io/pdbtopologyreader.cc \
    modules/io/gmxtrajectorywriter.cc  modules/io/pdbwriter.cc

libvotca_csg_la_LIBADD = $(VOTCA_TOOLS_LIBS) $(GMX_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)

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
