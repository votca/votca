.PHONY: programs scripts

all: cleantex programs.tex scripts.tex
	$(MAKE) $(MFLAGS) -C xml

%.tex:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $*
	ls -1 $*/*.tex | sed 's%\(.*\)\.tex%\\input{reference/\1}%' > $@

cleantex:
	rm -f programs.tex scripts.tex

%_subclean:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $* clean

clean: cleantex programs_subclean scripts_subclean
	$(MAKE) $(MFLAGS) -C xml clean
	rm -f *~
