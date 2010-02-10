AM_CPPFLAGS = $(CFLAGS_XML) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) -I$(top_srcdir)/include 

bin_PROGRAMS = \
    csg_map csg_dump csg_property csg_resample csg_stat csg_fmatch csg_imcrepack csg_gmxtopol # csg_nemat csg_conj csg_sample

csg_dump_SOURCES = csg_dump.cc

# csg_nemat_SOURCES = \
#    csg_nemat.cc

# csg_ndens_SOURCES = csg_ndens.cc

# csg_update_SOURCES = csg_update.cc

# csg_conj_SOURCES = csg_conj.cc

csg_resample_SOURCES = csg_resample.cc

csg_stat_SOURCES = csg_stat.cc imc.cc imc.h

csg_property_SOURCES = csg_property.cc

csg_fmatch_SOURCES = csg_fmatch.cc csg_fmatch_main.cc csg_fmatch.h

csg_imcrepack_SOURCES = csg_imcrepack.cc

csg_gmxtopol_SOURCES = csg_gmxtopol.cc

csg_map_SOURCES = csg_map.cc

# csg_sample_SOURCES = \
#    csg_sample.cc

AM_LDFLAGS =  

LDADD = ../libcsg/libvotca_csg.la $(LIBS_XML) ${LIBS_VOTCA_TOOLS} ${BOOST_LDFLAGS} ${BOOST_PROGRAM_OPTIONS_LIB} $(LIBS_GMX) ${GSL_LIBS} -lpthread
