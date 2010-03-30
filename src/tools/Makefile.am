AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(BOOST_CFLAGS)

bin_PROGRAMS = \
	ctp_map

ctp_map_SOURCES = csg_map.cc

#everything is already in libvotca_csg
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../libmd2qm/libvotca_md2qm.la
