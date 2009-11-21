AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

bin_PROGRAMS = \
    csg_boltzmann

csg_boltzmann_SOURCES = \
    bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc 

AM_LDFLAGS =  

LDADD = ../libcsg/libcsg.a  ${BOOST_PROGRAM_OPTIONS_LIB} $(LIBS_TOOLS) $(LIBS_GMX) $(LIBS_XML) $(LIBS_FFTW) -lm -lpthread 
