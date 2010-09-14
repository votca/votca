AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(VOTCA_MOO_CFLAGS) $(VOTCA_CSG_CFLAGS) $(VOTCA_KMC_CFLAGS) 

bin_PROGRAMS = \
	ctp_map ctp_printqmattop ctp_polymer ctp_run

ctp_map_SOURCES = ctp_map.cc md2qm_observer.cc
ctp_run_SOURCES = ctp_run.cc
ctp_printqmattop_SOURCES = printQMAtTop.cc atqmtopobserver.cc

#special treatment here since it's the only app which needs libkmc
#ctp_kmc_CPPFLAGS=$(VOTCA_KMC_CFLAGS) -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(VOTCA_MOO_CFLAGS) $(VOTCA_CSG_CFLAGS)
#ctp_kmc_LDFLAGS=$(VOTCA_KMC_LIBS) $(BIN_LDFLAGS)

ctp_polymer_SOURCES = polmain.cc polymer.cc
ctp_polymer_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(VOTCA_MOO_CFLAGS) $(VOTCA_CSG_CFLAGS) $(GSL_CFLAGS)
ctp_polymer_LDFLAGS = $(GSL_LIBS)

#everything is already in libvotca_csg
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../libmd2qm/libvotca_md2qm.la
