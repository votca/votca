
AM_CPPFLAGS = -I../../libcsg -I/usr/include/libxml2
AM_LDFLAGS = -L/misc/home/victorr/gmx/lib \
             -L../../libcsg/ \
             -lcsg \
             -lgmx -lboost_program_options -lxml2 -lm

bin_PROGRAMS = \
    csg_nemat csg_sample

csg_nemat_SOURCES = \
    csg_nemat.cc

csg_sample_SOURCES = \
    csg_sample.cc
