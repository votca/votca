AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

bin_PROGRAMS = \
    csg

csg_SOURCES = \
    bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc 

AM_LDFLAGS = ${BOOST_PROGRAM_OPTIONS_LIB}  

LDADD = ../libcsg/libcsg.a -ltools -lgmx -lxml2 -lfftw3 -lm
