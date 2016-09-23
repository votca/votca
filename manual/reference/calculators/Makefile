SHELL=/bin/bash

DEPS=$(wildcard $(XTP_PROGS))
calculators.tex: $(DEPS)
	./extract_calculators.sh $(XTP_PROGS)

dummy: ;

clean:                                 
	rm -f *~ calculators.tex *.backup

