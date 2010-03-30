AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(VOTCA_MOO_CFLAGS) $(VOTCA_CSG_CFLAGS) $(BOOST_CFLAGS)

bin_PROGRAMS = \
	ctp_map ctp_kmc

ctp_map_SOURCES = ctp_map.cc md2qm_observer.cc

ctp_kmc_SOURCES = kmc_cont.cc kmc_cont_app.cc
# special treatment here since it's the only app which needs libkmc
ctp_kmc_CPPFLAGS=$(VOTCA_KMC_CFLAGS) -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(VOTCA_MOO_CFLAGS) $(VOTCA_CSG_CFLAGS) $(BOOST_CFLAGS)
ctp_kmc_LDFLAGS=$(VOTCA_KMC_LIBS) $(BIN_LDFLAGS)

#everything is already in libvotca_csg
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../libmd2qm/libvotca_md2qm.la
