AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) $(BOOST_CFLAGS)

bin_PROGRAMS = \
   moo_overlap 

moo_overlap_SOURCES = moo_overlap.cc

AM_LDFLAGS = $(BIN_LDFLAGS) $(BOOST_LIBS)  $(VOTCA_TOOLS_LIBS)
LDADD = ../libmoo/libvotca_moo.la
