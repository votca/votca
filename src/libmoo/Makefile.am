AM_CPPFLAGS = -I$(top_srcdir)/include $(BOOST_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS)
lib_LTLIBRARIES = libvotca_moo.la

libvotca_moo_la_SOURCES = \
    basis_set.cpp  charges.cpp  crgunit.cc  crgunittype.cc \
    fock_matrix.cpp  mol_pair.cc  orbitals.cpp  qm_molecule.cpp

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_moo.pc

CLEANFILES = *~ 
