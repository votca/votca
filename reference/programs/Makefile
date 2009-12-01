DOCS=\
	csg_boltzmann.tex\
	csg_conj.tex\
	csg_fmatch.tex\
	csg_imcrepack.tex\
	csg_nemat.tex\
	csg_resample.tex\
	csg_call.tex\
	csg_dump.tex\
	csg_gmxtopol.tex\
	csg_inverse.tex\
	csg_property.tex\
	csg_stat.tex

SHELL=/bin/sh

.SUFFIXES: .t2t .tex
.PHONY: clean

all: $(DOCS)

%.t2t:
	@if [ -z "${CSGSHARE}" ]; then \
	  echo CSGSHARE not defined >&2; \
	  exit 1; \
	fi
	@if ! which $* &> /dev/null; then\
	  echo $* not found;\
	  exit 1;\
	fi
	@if ! $* --help &> /dev/null; then\
	  echo $* --help failed;\
	  exit 1;\
	fi
	echo -e "\n\n\n==$*==" > $@
	@#$* --help | sed 's/^  --/: --/;/--/s/[[:space:]]\{5,\}/&\n/' >> $@
	$* --help | sed 's/^[[:space:]]*-/- -/' >> $@

.t2t.tex:
	./txt2tags -t tex --no-header --infile $< --outfile $@

#cleant2t:
#	rm -f *.t2t

clean: 
	rm -f $(DOCS)
	rm -f *.tex *~
