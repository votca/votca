AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(CFLAGS_XML) -I$(top_srcdir)/include 

bin_PROGRAMS = \
    csg_dump csg_nemat csg_conj csg_property csg_resample csg_stat csg_fmatch csg_imcrepack # csg_sample

csg_dump_SOURCES = csg_dump.cc

csg_nemat_SOURCES = \
    csg_nemat.cc

# csg_ndens_SOURCES = csg_ndens.cc

# csg_update_SOURCES = csg_update.cc

csg_conj_SOURCES = csg_conj.cc

csg_resample_SOURCES = csg_resample.cc

csg_stat_SOURCES = csg_stat.cc imc.cc

csg_property_SOURCES = csg_property.cc

csg_fmatch_SOURCES = csg_fmatch.cc csg_fmatch_main.cc

csg_imcrepack_SOURCES = csg_imcrepack.cc

# csg_sample_SOURCES = \
#    csg_sample.cc

AM_LDFLAGS =  

LDADD = ../libcsg/libcsg.a -ltools ${BOOST_PROGRAM_OPTIONS_LIB} $(LIBS_GMX) $(LIBS_XML) -lgsl -lgslcblas -lm 
