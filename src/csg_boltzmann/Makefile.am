AM_CPPFLAGS = -I$(top_srcdir)/include

bin_PROGRAMS = \
	csg_boltzmann

csg_boltzmann_CPPFLAGS=$(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(BOOST_CFLAGS)

csg_boltzmann_SOURCES = \
	bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc \
	analysistool.h bondedstatistics.h stdanalysis.h tabulatedpotential.h

#everything is alread in libvotca_csg
csg_boltzmann_LDADD = ../libcsg/libvotca_csg.la
