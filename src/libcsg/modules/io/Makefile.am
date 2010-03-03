AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_LTLIBRARIES = libio.la

libio_la_SOURCES = \
	gmxtopologyreader.cc 	gmxtopologyreader.h 	gmxtrajectoryreader.cc\
	gmxtrajectoryreader.h 	gmxtrajectorywriter.cc 	gmxtrajectorywriter.h\
	grotopologyreader.cc 	grotopologyreader.h 	growriter.cc\
	growriter.h      	\
	pdbtopologyreader.cc 	pdbtopologyreader.h 	pdbwriter.cc\
	pdbwriter.h 		xmltopologyreader.cc	xmltopologyreader.h

libio_la_LIBADD = $(GMX_LIBS) $(VOTCA_TOOLS_LIBS)
libio_la_CPPFLAGS = $(GMX_CFLAGS) $(VOTCA_TOOLS_CFLAGS) $(AM_CPPFLAGS)

CLEANFILES = *~
