.PHONY: programs

all: programs.tex
	$(MAKE) $(MFLAGS) -C src

programs.tex: programs 
	ls -1 programs/*.tex | sed 's%\(.*\)\.tex%\\input{reference/\1}%' > programs.tex

programs:
	$(MAKE) $(MFLAGS) -C programs

clean:
	rm -f programs.tex
	$(MAKE) $(MFLAGS) -C src clean
	$(MAKE) $(MFLAGS) -C programs clean
