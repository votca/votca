AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

lib_LIBRARIES = libcsg.a

libcsg_a_LIBADD = modules/io/libcsgio.a
libcsg_a_DEPENDENCIES = modules/io/libcsgio.a

libcsg_a_SOURCES = \
    cgengine.cc       exclusionlist.cc    neighbourlist.cc  topologymap.cc \
    cgmoleculedef.cc  imccoefficients.cc  nematicorder.cc   topologyreader.cc \
    configuration.cc  map.cc              numberdist.cc     trajectoryreader.cc \
    connectivity.cc   moleculeinfo.cc     topology.cc       trajectorywriter.cc \
    modules/io/gmxtopologyreader.cc    modules/io/grotopologyreader.cc  modules/io/mdptrajectoryreader.cc  modules/io/xmltopologyreader.cc \
    modules/io/gmxtrajectoryreader.cc  modules/io/growriter.cc          modules/io/pdbtopologyreader.cc \
    modules/io/gmxtrajectorywriter.cc  modules/io/mdptopologyreader.cc  modules/io/pdbwriter.cc

CLEANFILES = *~
#DISTCLEANFILES =
