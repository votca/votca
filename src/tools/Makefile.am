AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(BOOST_CFLAGS)

bin_PROGRAMS = \
	ctp_map ctp_integrals ctp_ratecalc ctp_cont_kmc

ctp_map_SOURCES = atqmtopobserver.cc ctp_map.cc intergalsapp.cc integrals.cc kmc_cont_app.cpp kmc_cont.cpp md2qm_observer.cc printQMAtTop.cc proj.cc projobserver.cc qmapplication.cc ratecalcapp.cc ratecalc.cc

#everything is already in libvotca_csg
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../libmd2qm/libvotca_md2qm.la
