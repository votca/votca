PICS=\
	flowchart.eps\

all: $(PICS)

.SUFFIXES: .tex .eps
.PHONY: clean

.tex.eps:
	./tex2eps $< $@

clean:
	rm -f $(PICS)
	rm -f *~
