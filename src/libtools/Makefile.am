AM_CPPFLAGS = -I$(top_srcdir)/include
lib_LTLIBRARIES = libvotca_tools.la

BUILT_SOURCES = hgversion.h

libvotca_tools_la_LDFLAGS = -no-undefined $(GSL_LIBS) $(FFTW_LIBS) $(EXPAT_LIBS) $(BOOST_LIBS) $(PTHREAD_LIBS) -lm
libvotca_tools_la_CPPFLAGS = $(AM_CPPFLAGS) $(GSL_CFLAGS) $(FFTW_CFLAGS) $(EXPAT_CFLAGS) $(BOOST_CFLAGS) $(PTHREAD_CFLAGS)

libvotca_tools_la_SOURCES = \
application.cc correlate.cc       datacollection.cc  matrix.cc    random.cc       tokenizer.cc \
crosscorrelate.cc  histogram.cc       parcer.cc    rangeparser.cc \
cubicspline.cc     histogramnew.cc    property.cc  table.cc version.cc linalg.cc parsexml.cc \
mutex.cc thread.cc

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
dummy: ;
