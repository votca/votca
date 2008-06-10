
AM_CPPFLAGS = -I../../libcsg -I/usr/include/libxml2

bin_PROGRAMS = \
    csg

csg_SOURCES = \
    main.cc  stdanalysis.cc  tabulatedpotential.cc 

AM_LDFLAGS = -L/misc/home/victorr/gmx/lib \
             -L../../libcsg/dist 

LDADD = -lcsg -lgmx  -lboost_program_options -lxml2 -lfftw3 -lm
