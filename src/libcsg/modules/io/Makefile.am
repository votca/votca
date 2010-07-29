AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_LTLIBRARIES = libio.la

if GMX
GMX_SOURCES = \
	gmxtopologyreader.cc 	gmxtopologyreader.h 	gmxtrajectoryreader.cc\
	gmxtrajectoryreader.h 	gmxtrajectorywriter.cc 	gmxtrajectorywriter.h\
	grotopologyreader.cc 	grotopologyreader.h 	growriter.cc\
	growriter.h      	pdbtopologyreader.cc 	pdbtopologyreader.h
endif

libio_la_SOURCES = $(GMX_SOURCES)\
	esptrajectoryreader.cc esptrajectoryreader.h \
	esptopologyreader.cc    esptopologyreader.h  	pdbwriter.cc\
	pdbwriter.h 		xmltopologyreader.cc	xmltopologyreader.h

libio_la_LDFLAGS = -no-undefined
libio_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(GMX_CFLAGS)

CLEANFILES = *~
