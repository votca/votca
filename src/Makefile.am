
AM_CPPFLAGS = -I../../libcsg -I/usr/include/libxml2

bin_PROGRAMS = \
    csg_dump csg_nemat csg_ndens csg_update # csg_sample

csg_dump_SOURCES = csg_dump.cc

csg_nemat_SOURCES = \
    csg_nemat.cc

csg_ndens_SOURCES = csg_ndens.cc

csg_update_SOURCES = csg_update.cc

# csg_sample_SOURCES = \
#    csg_sample.cc

AM_LDFLAGS = -L/misc/home/victorr/gmx/lib \
             -L../../libcsg/dist 

LDADD = -lcsg -lgmx  -lboost_program_options -lxml2 -lm
