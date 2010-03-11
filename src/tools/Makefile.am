AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(BOOST_CFLAGS)

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

#fmatch uses gsl
csg_fmatch_SOURCES = csg_fmatch.cc csg_fmatch_main.cc csg_fmatch.h
csg_fmatch_CPPFLAGS = $(AM_CPPFLAGS) $(GSL_CFLAGS)
csg_fmatch_LDFLAGS = $(AM_LDFLAGS)

csg_imcrepack_SOURCES = csg_imcrepack.cc

csg_gmxtopol_SOURCES = csg_gmxtopol.cc

csg_map_SOURCES = csg_map.cc

# csg_sample_SOURCES = \
#    csg_sample.cc

#everything is already in libvotca_csg
AM_LDFLAGS = $(BIN_LDFLAGS) ../libcsg/libvotca_csg.la
