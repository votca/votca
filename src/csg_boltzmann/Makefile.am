AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

bin_PROGRAMS = \
	csg_boltzmann

csg_boltzmann_SOURCES = \
	bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc \
	analysistool.h bondedstatistics.h stdanalysis.h tabulatedpotential.h

AM_LDFLAGS =  

LDADD = ../libcsg/libvotca_csg.a  ${BOOST_PROGRAM_OPTIONS_LIB} $(LIBS_VOTCA_TOOLS) $(LIBS_GMX) $(LIBS_XML) $(LIBS_FFTW) -lm -lpthread 
