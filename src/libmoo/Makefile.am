@SET_MAKE@
ACLOCAL_AMFLAGS = -I config
AUTOMAKE_OPTIONS = foreign

AM_LDFLAGS = $(VOTCA_TOOLS_LIBS)

AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS) 
lib_LTLIBRARIES = libvotca_moo.la

libvotca_moo_la_SOURCES = \
    basis_set.cpp  charges.cpp  crgunit.cc  crgunittype.cc \
    fock.cc  orbitals.cpp  mol_and_orb.cpp jcalc.cc
    
CLEANFILES = *~ *.a
