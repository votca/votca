HGVERSION:= $(shell hg -R $(top_srcdir) parents --template 'hgid: {node|short}' 2> /dev/null || grep node $(top_srcdir)/.hg_archival.txt 2> /dev/null || true )
AM_CPPFLAGS = $(CFLAGS_XML) $(GSL_CFLAGS) -I$(top_srcdir)/include -DHGVERSION="\"${HGVERSION}\""
lib_LTLIBRARIES = libvotca_tools.la

libvotca_tools_la_LIBADD = $(LIBS_XML) $(GSL_LIBS)

libvotca_tools_la_SOURCES = \
correlate.cc       datacollection.cc  matrix.cc    random.cc       tokenizer.cc \
crosscorrelate.cc  histogram.cc       parcer.cc    rangeparser.cc \
cubicspline.cc     histogramnew.cc    property.cc  table.cc version.cc

CLEANFILES = *~ *.a hgstamp

version.o: hgstamp
version.lo: hgstamp
hgstamp: dummy 
	[ -f $@ ] || touch $@
	echo $(HGVERSION) | cmp -s $@ - || echo $(HGVERSION) > $@
dummy: ;
