templatedir = ${datadir}/votca/template

template_DATA = template.cc Makefile.template CMakeLists.txt.template

template.$(OBJEXT): votca/csg

votca/csg:
	$(MKDIR_P) votca
	cd votca && rm -f csg && $(LN_S) $(abs_top_builddir)/include csg

install-data-hook:
	cd $(DESTDIR)$(templatedir) && mv Makefile.template Makefile
	cd $(DESTDIR)$(templatedir) && mv CMakeLists.txt.template CMakeLists.txt

clean-local:
	rm -rf votca

AM_CPPFLAGS = $(VOTCA_TOOLS_CFLAGS)
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../../src/libcsg/libvotca_csg.la

template_SOURCES = template.cc

noinst_PROGRAMS = template

CLEANFILES = *~

