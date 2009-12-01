.PHONY: programs bootstrap_scripts

all: bootstrap_scripts cleantex programs.tex scripts.tex
	$(MAKE) $(MFLAGS) -C src

bootstrap_scripts:
	$(MAKE) $(MFLAGS) -C scripts -f Makefile.in bootstrap

%.tex:
	$(MAKE) $(MFLAGS) -C $*
	ls -1 $*/*.tex | sed 's%\(.*\)\.tex%\\input{reference/\1}%' > $@

cleantex:
	rm -f programs.tex scripts.tex

clean: cleantex bootstrap_scripts
	$(MAKE) $(MFLAGS) -C src clean
	$(MAKE) $(MFLAGS) -C programs clean
	$(MAKE) $(MFLAGS) -C scripts clean
	rm -f *~
