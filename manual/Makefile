SHELL=/bin/bash
VER=1.3-dev
HGID:=$(shell hg parents -R . --template "{node|short}" | sed 's/.*/\\newcommand{\\hgid}{${VER} (&)}/')
LATEXMK=./scripts/latexmk.pl
LATEXMKOPTS=-e '$$latex=q/latex --halt-on-error %O %S/'

NAME=manual
ifeq ($(OSTYPE),darwin)
CWD=$(shell pwd)
ND=$(subst work/votca_ct,Dropbox,$(CWD))
endif

all: $(NAME).pdf
dvi: $(NAME).dvi
ps: $(NAME).ps

$(NAME).tex: hgid.tex reference_submake fig_submake programs_submake theory_submake reference_submake input_submake appendix_submake titlepage.tex manual.bib           

#remove broken dvi if LATEXMK fails
.DELETE_ON_ERROR: %.dvi

%.dvi: %.tex dummy
	./abbreviate_journals.pl 
	$(LATEXMK) $(LATEXMKOPTS) -dvi $<

%.pdf: %.dvi
ifeq ($(OSTYPE),darwin)
	dvips $(NAME).dvi
	ps2pdf $(NAME).ps
	mkdir -p $(ND)
	cp $(NAME).pdf $(ND)
else
	dvipdf $*
endif


%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	rm -f *.backup
	$(MAKE) $(MFLAGS) -C $* clean

qclean:
	$(LATEXMK) -C $(NAME).tex

clean: qclean reference_subclean fig_subclean programs_subclean theory_subclean reference_subclean input_subclean
	rm -f $(NAME).fdb_latexmk $(NAME).brf literature_short.bib literature.bib
	rm -f hgid.tex
	rm -f *~

tar: all
	rm -f $(NAME).tar.gz
	#dirty sed hack ahead to grep files used by latexmk to produce the pdf
	tar cvzhf $(NAME).tar.gz $(NAME).pdf \
  	   `sed -n 's/^[[:space:]]*"\([^/][^"]*\.\(bib\|tex\|eps\|cls\)\)".*$$/\1/p' $(NAME).fdb_latexmk` 

hgid.tex: dummy
	[ -f hgid.tex ] || touch hgid.tex
	echo '$(HGID)' | cmp -s hgid.tex - || echo '$(HGID)' > hgid.tex

dummy: ;
