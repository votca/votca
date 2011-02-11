
templatedir = ${datadir}/votca/template

template_DATA = template.cc template_threaded.cc Makefile.template CMakeLists.txt.template

BUILT_SOURCES = votca/csg

votca/csg:
	mkdir votca
#this will work for global and local paths in top_srcdir
	incl=`$(am__cd) $(top_srcdir)/include; pwd`; \
	cd votca && $(LN_S) $$incl csg

install-data-hook:
	cd $(DESTDIR)$(templatedir) && mv Makefile.template Makefile
	cd $(DESTDIR)$(templatedir) && mv CMakeLists.txt.template CMakeLists.txt

uninstall-hook:
	cd $(DESTDIR)$(templatedir) && rm -f Makefile
	cd $(DESTDIR)$(templatedir) && rm -f CMakeLists.txt

clean-local:
	rm -rf votca

AM_CPPFLAGS = $(VOTCA_TOOLS_CFLAGS)
AM_LDFLAGS = $(BIN_LDFLAGS)
LDADD = ../../src/libcsg/libvotca_csg.la

template_SOURCES = template.cc
template_threaded_SOURCES = template_threaded.cc

noinst_PROGRAMS = template template_threaded

EXTRA_DIST = $(template_DATA)

CLEANFILES = *~

