AM_CPPFLAGS = -I$(top_srcdir)/include

noinst_LTLIBRARIES = libio.la

if GMX
GMX_SOURCES = gmx_libs_version.h \
	gmxtopologyreader.cc 	gmxtopologyreader.h 	gmxtrajectoryreader.cc\
	gmxtrajectoryreader.h 	gmxtrajectorywriter.cc 	gmxtrajectorywriter.h\
	grotopologyreader.cc 	grotopologyreader.h 	growriter.cc\
	growriter.h      	pdbtopologyreader.cc 	pdbtopologyreader.h \
	gmx_version.cc gmx_version.h

noinst_PROGRAMS = gmx_print_version

gmx_print_version_SOURCES = gmx_print_version.cc
gmx_print_version_CPPFLAGS = $(GMX_CFLAGS)
gmx_print_version_LDFLAGS= $(GMX_LIBS)

BUILT_SOURCES = gmx_libs_version.h

gmx_libs_version.h: gmx_print_version
	./$< > $@ || echo 'static const std::string gmx_version = "UNDEFINED";' > $@
	
endif

libio_la_SOURCES = $(GMX_SOURCES)\
	esptrajectoryreader.cc esptrajectoryreader.h \
	esptopologyreader.cc    esptopologyreader.h  	pdbwriter.cc\
	pdbwriter.h 		xmltopologyreader.cc	xmltopologyreader.h \
	lammpsreader.h lammpsreader.cc

libio_la_LDFLAGS = -no-undefined
libio_la_CPPFLAGS = $(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS) $(GMX_CFLAGS)

CLEANFILES = *~
