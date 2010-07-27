templatedir = ${datadir}/votca/template

template_DATA = template.cc Makefile.template 

install-data-hook:
	cd $(DESTDIR)$(templatedir) && mv Makefile.template Makefile

AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS)
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../../src/libcsg/libvotca_csg.la

template_SOURCES = template.cc

#noinst_PROGRAMS = template

CLEANFILES = *~

