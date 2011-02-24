AM_CPPFLAGS = -I$(top_srcdir)/include

bin_PROGRAMS = \
	csg_boltzmann

if HAVE_TXT2TAGS
man1_MANS= csg_boltzmann.man
%.man: %.t2t
	${TXT2TAGS} -t man -i $< -o $@

%.t2t: %.help $(top_builddir)/scripts/help2t2t
	$(top_builddir)/scripts/help2t2t $< > $@

%.help: %
	[ -z "$${LD_LIBRARY_PATH}" ] && export LD_LIBRARY_PATH="../libcsg/.libs" || export LD_LIBRARY_PATH="../libcsg/.libs:$$LD_LIBRARY_PATH"; \
	./$* --help > $@

else
man1_MANS=
endif
csg_boltzmann_CPPFLAGS=$(AM_CPPFLAGS) $(VOTCA_TOOLS_CFLAGS)

csg_boltzmann_SOURCES = \
analysistool.h bondedstatistics.cc bondedstatistics.h main.cc \
stdanalysis.cc stdanalysis.h tabulatedpotential.cc tabulatedpotential.h

#everything is alread in libvotca_csg
csg_boltzmann_LDADD = ../libcsg/libvotca_csg.la
csg_boltzmann_LDFLAGS = $(BIN_LDFLAGS)

CLEANFILES = $(man1_MANS) *~
