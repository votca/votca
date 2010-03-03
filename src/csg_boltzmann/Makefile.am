bin_PROGRAMS = \
	csg_boltzmann

csg_boltzmann_CPPFLAGS=$(BOOST_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) -I$(top_srcdir)/include

csg_boltzmann_SOURCES = \
	bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc \
	analysistool.h bondedstatistics.h stdanalysis.h tabulatedpotential.h

csg_boltzmann_LDADD = ../libcsg/libvotca_csg.la  ${BOOST_LDFLAGS} ${BOOST_PROGRAM_OPTIONS_LIB} $(LIBS_VOTCA_TOOLS) -lm
