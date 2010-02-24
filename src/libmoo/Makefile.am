@SET_MAKE@
ACLOCAL_AMFLAGS = -I config
AUTOMAKE_OPTIONS = foreign


AM_CPPFLAGS = -I$(top_srcdir)/include -I../../../libcsg 
lib_LTLIBRARIES = libvotca_moo.la

libvotca_moo_la_SOURCES = \
    basis_set.cpp  charges.cpp  crgunit.cc  crgunittype.cc \
    fock_matrix.cpp  mol_pair.cc  orbitals.cpp  qm_molecule.cpp
    
CLEANFILES = *~ *.a
