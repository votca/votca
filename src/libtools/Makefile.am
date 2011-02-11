AM_CPPFLAGS = -I$(top_srcdir)/include
lib_LTLIBRARIES = libvotca_tools.la

BUILT_SOURCES = hgversion.h

libvotca_tools_la_LDFLAGS = -no-undefined $(GSL_LIBS) $(FFTW_LIBS) $(EXPAT_LIBS) $(BOOST_LIBS) $(PTHREAD_LIBS) -lm
libvotca_tools_la_CPPFLAGS = $(AM_CPPFLAGS) $(GSL_CFLAGS) $(FFTW_CFLAGS) $(EXPAT_CFLAGS) $(BOOST_CFLAGS) $(PTHREAD_CFLAGS)

libvotca_tools_la_SOURCES = \
akimaspline.cc application.cc correlate.cc crosscorrelate.cc \
cubicspline.cc datacollection.cc easylock.cc histogram.cc histogramnew.cc \
linalg.cc linspline.cc matrix.cc mutex.cc parcer.cc parsexml.cc \
property.cc random.cc rangeparser.cc spline.cc table.cc thread.cc \
tokenizer.cc version.cc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_tools.pc

CLEANFILES = *~ *.a hgversion.h

install-exec-hook:
if NO_LA_FILES
	rm -f $(DESTDIR)$(libdir)/libvotca_tools.la
endif

HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_srcdir)/.hg_archival.txt 2> /dev/null || true )
hgversion.h: dummy
	[ -f $@ ] || touch $@
	echo 'static const std::string hgversion = "$(HGVERSION)";' | cmp -s $@ - \
	  || echo 'static const std::string hgversion = "$(HGVERSION)";' > $@

.PHONY: dummy
dummy: ;
