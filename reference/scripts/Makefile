SCRIPTS=\
	functions_inverse.sh.tex\
	initalize_generic.sh.tex\
	initalize_generic.sh.tex\
	update_ibm.sh.tex\
	update_imc.sh.tex\
	add_pot_ibm.sh.tex\
	add_pot_ibm.sh.tex\
	RDF_to_POT.pl.tex\
	post_update.sh.tex\
	post_update_single.sh.tex\
	postupd_pressure.sh.tex\
	postupd_splinesmooth.sh.tex\
	postupd_smooth.sh.tex\
	post_add.sh.tex\
	post_add_single.sh.tex\
	dpot_shift_nb.pl.tex\
	resample_to_calc.sh.tex\
	dpot_crop.pl.tex\
	update_ibm_single.sh.tex\
	update_ibm_pot.pl.tex\
	update_imc_single.sh.tex\
	solve_matlab.sh.tex\
	solve_octave.sh.tex\
	imc_purify.sh.tex\
	pressure_cor_simple.pl.tex\
	pressure_cor_wjk.pl.tex\
	add_POT.pl.tex\
	integrate.pl.tex\
	extrapolate.pl.tex\
	merge_tables.pl.tex\
	smooth.pl.tex\
	table_linearop.pl.tex\
	tables_jackknife.pl.tex\
	initalize_gromacs.sh.tex\
	prepare_gromacs.sh.tex\
	run_gromacs.sh.tex\
	calc_pressure_gromacs.sh.tex\
	calc_rdf_gromacs.sh.tex\
	imc_stat_generic.sh.tex\
	potential_to_gromacs.sh.tex\
	table_to_xvg.pl.tex\
	functions_gromacs.sh.tex\

.SUFFIXES: .t2t .tex
.PHONY: clean

SHELL=/bin/sh

all: $(SCRIPTS)

%.t2t:
	@if [ -z "${CSGSHARE}" ]; then \
	  echo CSGSHARE not defined >&2; \
	  exit 1; \
	fi
	@if ! cd ${CSGSHARE}/scripts/inverse && ./$* --help &> /dev/null; then\
	  echo $* --help failed;\
	  exit 1;\
	fi
	echo -e "\n\n\n==$*==" > $@
	@#$* --help | sed 's/^  --/: --/;/--/s/[[:space:]]\{5,\}/&\n/' >> $@
	( cd ${CSGSHARE}/scripts/inverse && ./$* --help ) | sed 's/^[[:space:]]*-/- -/' >> $@

.t2t.tex:
	./txt2tags -t tex --no-header --infile $< --outfile $@

#cleant2t:
#	rm -f *.t2t

clean: 
	rm -f $(SCRIPTS)

Makefile: Makefile.in
	@if [ -z "${CSGSHARE}" ]; then \
	  echo CSGSHARE not defined >&2; \
	  exit 1; \
	fi
	cat ${CSGSHARE}/scripts/inverse/csg_table | awk '/^[^#]/{print $$3}' | \
	  sed -e 's/.*/\t&.tex\\/' -e '1i SCRIPTS=\\' > $@
	echo >> $@
	cat Makefile.in >> $@
	
