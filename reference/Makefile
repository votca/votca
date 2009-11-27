.PHONY: programs

all: programs.tex scripts.tex
	$(MAKE) $(MFLAGS) -C src

%.tex:
	$(MAKE) $(MFLAGS) -C $*
	ls -1 $*/*.tex | sed 's%\(.*\)\.tex%\\input{reference/\1}%' > $@

clean:
	rm -f programs.tex scripts.tex
	$(MAKE) $(MFLAGS) -C src clean
	$(MAKE) $(MFLAGS) -C programs clean
	$(MAKE) $(MFLAGS) -C scripts clean
