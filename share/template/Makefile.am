templatedir = ${datadir}/votca/template

template_DATA = csg_template.cc Makefile.template 

install-data-hook:
	cd $(DESTDIR)$(templatedir) && mv Makefile.template Makefile

AM_CPPFLAGS = -I$(top_srcdir)/include $(VOTCA_TOOLS_CFLAGS)
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../../src/libcsg/libvotca_csg.la

csg_template_SOURCES = csg_template.cc

noinst_PROGRAMS = csg_template

CLEANFILES = *~

