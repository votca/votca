AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

bin_PROGRAMS = \
    csg

csg_SOURCES = \
    bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc 

AM_LDFLAGS =  

LDADD = ../libcsg/libcsg.a  ${BOOST_PROGRAM_OPTIONS_LIB} -ltools $(LIBS_GMX) -lxml2 -lfftw3 -lpthread -lz -lm
