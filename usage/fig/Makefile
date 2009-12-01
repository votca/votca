PICS=\
	flow_boltzmann.eps\
	flow_fmatch.eps\
	flow_ibi.eps\

all: $(PICS)

.SUFFIXES: .tex .eps
.PHONY: clean

.tex.eps:
	./tex2eps $< $@

clean:
	rm -f $(PICS)
	rm -f *~
