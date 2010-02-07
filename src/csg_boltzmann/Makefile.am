AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

bin_PROGRAMS = \
	csg_boltzmann

csg_boltzmann_SOURCES = \
	bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc \
	analysistool.h bondedstatistics.h stdanalysis.h tabulatedpotential.h

csg_boltzmann_LDADD = ../libcsg/libvotca_csg.la  ${BOOST_LDFLAGS} ${BOOST_PROGRAM_OPTIONS_LIB} $(LIBS_VOTCA_TOOLS) $(LIBS_GMX) $(GSL_LIBS) $(LIBS_XML) $(LIBS_FFTW) -lm -lpthread
