AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

bin_PROGRAMS = \
    csg

csg_SOURCES = \
    bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc 

AM_LDFLAGS =  

LDADD = ../libcsg/libcsg.a  ${BOOST_PROGRAM_OPTIONS_LIB} -ltools -lgmx_mpi -lmpi -llam -lxml2 -lfftw3 -lm
