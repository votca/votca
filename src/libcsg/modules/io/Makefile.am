AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_LTLIBRARIES = libio.la

libio_la_SOURCES = \
	gmxtopologyreader.cc 	gmxtopologyreader.h 	gmxtrajectoryreader.cc\
	gmxtrajectoryreader.h 	gmxtrajectorywriter.cc 	gmxtrajectorywriter.h\
	grotopologyreader.cc 	grotopologyreader.h 	growriter.cc\
	growriter.h      	\
	pdbtopologyreader.cc 	pdbtopologyreader.h 	pdbwriter.cc\
	pdbwriter.h 		xmltopologyreader.cc	xmltopologyreader.h

#no LIBADD here, because this is a noint lib
libio_la_LIBADD =
libio_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(BOOST_CFLAGS) $(GMX_CFLAGS) $(XML_CFLAGS)

CLEANFILES = *~
