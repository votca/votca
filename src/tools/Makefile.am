AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include 

bin_PROGRAMS = \
    csg_dump csg_nemat csg_ndens csg_update csg_conj # csg_fmatch # csg_sample

csg_dump_SOURCES = csg_dump.cc

csg_nemat_SOURCES = \
    csg_nemat.cc

csg_ndens_SOURCES = csg_ndens.cc

csg_update_SOURCES = csg_update.cc

csg_conj_SOURCES = csg_conj.cc

# csg_fmatch_SOURCES = csg_fmatch.cc

# csg_sample_SOURCES = \
#    csg_sample.cc

AM_LDFLAGS =  

LDADD = ../libcsg/libcsg.a -ltools ${BOOST_PROGRAM_OPTIONS_LIB} $(LIBS_GMX) -lxml2 -lgsl -lgslcblas -lm 
