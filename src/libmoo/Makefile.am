@SET_MAKE@
ACLOCAL_AMFLAGS = -I config
AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS)
lib_LTLIBRARIES = libvotca_moo.la

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_moo.pc

libvotca_moo_la_SOURCES = \
    basis_set.cpp  charges.cpp  crgunit.cc  crgunittype.cc \
    fock.cc  orbitals.cpp  mol_and_orb.cpp jcalc.cc
    
libvotca_moo_la_LDFLAGS = -no-undefined $(VOTCA_TOOLS_LIBS)

CLEANFILES = *~ *.a
