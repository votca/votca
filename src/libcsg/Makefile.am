HGVERSION:= $(shell grep -q node ../../hg_archival.txt 2> /dev/null || hg parents --template 'hgid: {node}' 2> /dev/null || true )
AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

lib_LIBRARIES = libvotca_csg.a

libvotca_csg_a_SOURCES = \
    beadlist.cc       exclusionlist.cc    nblist.cc         topology.cc          trajectorywriter.cc \
   cgengine.cc       imccoefficients.cc  neighbourlist.cc  topologymap.cc \
   cgmoleculedef.cc  map.cc              nematicorder.cc   topologyreader.cc \
   connectivity.cc   molecule.cc         numberdist.cc     trajectoryreader.cc \
   imcio.cc nblistgrid.cc version.cc \
   modules/io/gmxtopologyreader.cc    modules/io/grotopologyreader.cc  modules/io/xmltopologyreader.cc \
    modules/io/gmxtrajectoryreader.cc  modules/io/growriter.cc          modules/io/pdbtopologyreader.cc \
    modules/io/gmxtrajectorywriter.cc  modules/io/pdbwriter.cc

CLEANFILES = *~
#DISTCLEANFILES =
