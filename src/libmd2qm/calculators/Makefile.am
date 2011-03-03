AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_LTLIBRARIES = libcalculators.la

libcalculators_la_SOURCES = \
    contkmc.h  estatics.cc  estatics.h  integrals.h  marcusrates.h  readxml.cc  readxml.h  \
writexml.h shufflenrg.cc  shufflenrg.h generate_nrgs.cc  generate_nrgs.h energycorr.cc  energycorr.h \
polymerrates.h polymerrates.cc CurrentDensity.h

libcalculators_la_LDFLAGS = -no-undefined
libcalculators_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(VOTCA_CSG_CFLAGS) $(VOTCA_MOO_CFLAGS)

CLEANFILES = *~
