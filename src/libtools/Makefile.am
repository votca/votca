HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_srcdir)/.hg_archival.txt 2> /dev/null || true )
AM_CPPFLAGS = -I$(top_srcdir)/include -DHGVERSION="\"$(HGVERSION)\""
lib_LTLIBRARIES = libvotca_tools.la

#no BOOST_LDFLAGS here, because we just need the headers
libvotca_tools_la_LDFLAGS = -no-undefined $(GSL_LIBS) $(FFTW_LIBS) $(EXPAT_LIBS) -lm
libvotca_tools_la_CPPFLAGS = $(AM_CPPFLAGS) $(GSL_CFLAGS) $(FFTW_CFLAGS) $(EXPAT_CFLAGS) $(BOOST_CFLAGS)

libvotca_tools_la_SOURCES = \
correlate.cc       datacollection.cc  matrix.cc    random.cc       tokenizer.cc \
crosscorrelate.cc  histogram.cc       parcer.cc    rangeparser.cc \
cubicspline.cc     histogramnew.cc    property.cc  table.cc version.cc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libvotca_tools.pc

CLEANFILES = *~ *.a hgstamp

if NO_LA_FILES
install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/libvotca_tools.la
endif

version.o: hgstamp
version.lo: hgstamp
hgstamp: dummy 
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@
dummy: ;
