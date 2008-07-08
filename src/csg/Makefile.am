AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include

bin_PROGRAMS = \
    csg

csg_SOURCES = \
    bondedstatistics.cc main.cc  stdanalysis.cc  tabulatedpotential.cc 

AM_LDFLAGS =  

LDADD = ../libcsg/libcsg.a -ltools -lgmx  -lboost_program_options -lxml2 -lfftw3 -lm
