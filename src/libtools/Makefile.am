HGVERSION:= $(shell grep -q node ../../hg_archival.txt 2> /dev/null || hg parents --template 'hgid: {node}' 2> /dev/null || true )
AM_CPPFLAGS = $(CFLAGS_XML) -I$(top_srcdir)/include -DHGVERSION="\"${HGVERSION}\""
lib_LIBRARIES = libvotca_tools.a

libvotca_tools_a_SOURCES = \
correlate.cc       datacollection.cc  matrix.cc    random.cc       tokenizer.cc \
crosscorrelate.cc  histogram.cc       parcer.cc    rangeparser.cc \
cubicspline.cc     histogramnew.cc    property.cc  table.cc version.cc

CLEANFILES = *~ *.a
