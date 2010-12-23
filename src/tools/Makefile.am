AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS)

bin_PROGRAMS = \
    csg_map csg_dump csg_property csg_resample csg_stat csg_fmatch csg_imcrepack csg_gmxtopol csg_part_dist # csg_nemat csg_conj csg_sample

if HAVE_TXT2TAGS
man1_MANS= \
	csg_map.man csg_dump.man csg_property.man csg_resample.man csg_stat.man csg_fmatch.man\
	csg_imcrepack.man csg_gmxtopol.man csg_part_dist.man
%.man: %.t2t
	${TXT2TAGS} -t man -i $< -o $@

%.t2t: % $(top_builddir)/scripts/help2t2t
	$(top_builddir)/scripts/help2t2t $* > $@

else
man1_MANS=
endif

csg_dump_SOURCES = csg_dump.cc

# csg_nemat_SOURCES = \
#    csg_nemat.cc

# csg_ndens_SOURCES = csg_ndens.cc

# csg_update_SOURCES = csg_update.cc

# csg_conj_SOURCES = csg_conj.cc

csg_resample_SOURCES = csg_resample.cc

csg_stat_SOURCES = csg_stat.cc imc.cc imc.h

csg_property_SOURCES = csg_property.cc

csg_fmatch_SOURCES = csg_fmatch.cc csg_fmatch.h

csg_imcrepack_SOURCES = csg_imcrepack.cc

csg_gmxtopol_SOURCES = csg_gmxtopol.cc

csg_map_SOURCES = csg_map.cc

csg_part_dist_SOURCES = csg_part_dist.cc

# csg_sample_SOURCES = \
#    csg_sample.cc

#everything is already in libvotca_csg
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../libcsg/libvotca_csg.la

CLEANFILES = $(man1_MANS) *~

