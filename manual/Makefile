SHELL=/bin/bash

GITID_PLAIN:=$(shell git rev-parse --short HEAD 2> /dev/null)
GITID:=$(shell echo $(GITID_PLAIN) | sed -e 's/.*/\\newcommand{\\gitid}{${VER} (&)}/' -e 's/_rc/\\&/' -e 's/()//' )
LATEXMK=./scripts/latexmk.pl
LATEXMKOPTS=-halt-on-error --extra-mem-bot=100000000

NAME=xtp-manual
ifeq ($(OSTYPE),darwin)
CWD=$(shell pwd)
ND=$(subst work/votca_xtp,Dropbox,$(CWD))
endif

all: check $(NAME).pdf
dvi: $(NAME).dvi
ps: $(NAME).ps

$(NAME).tex: gitid.tex fig_submake programs_submake theory_submake reference_submake input_submake appendix_submake titlepage.tex $(NAME).bib

%.dvi: %.tex dummy
	./abbreviate_journals.pl 
	$(LATEXMK) $(LATEXMKOPTS) -dvi  $< || rm -vf $*.dvi
	@[ -f $@ ]

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

check:
	@if [ -z '$(XTP_SCRIPTS)' ] || [ -z '$(XTP_PROGS)' ]; then \
	  echo '#################################';\
	  echo '#                               #';\
	  echo '#  XTP_SCRIPTS or XTP_PROGS not #';\
	  echo '#  defined, call "make manual"  #';\
	  echo '#  in the main directory !!     #';\
	  echo '#                               #';\
	  echo '#################################';\
	  exit 1;\
	fi

install_check:
	@if [ -z '$(XTP_MANUAL_LOCATION)' ]; then \
	  echo '######################################';\
	  echo '#                                    #';\
	  echo '#  XTP_MANUAL_LOCATION not defined   #';\
	  echo '#  call "make installmanual" in the  #';\
	  echo '#  main directory !!                 #';\
	  echo '#                                    #';\
	  echo '######################################';\
	  exit 1;\
	fi
install: install_check $(NAME).pdf
	mkdir -p $(DESTDIR)$(XTP_MANUAL_LOCATION)
	echo "Installing to $(DESTDIR)${XTP_MANUAL_LOCATION}"
	cp -v $(NAME).pdf $(DESTDIR)$(XTP_MANUAL_LOCATION)/$(NAME).pdf

clean: qclean reference_subclean fig_subclean programs_subclean theory_subclean reference_subclean input_subclean
	rm -f $(NAME).fdb_latexmk $(NAME).brf literature_short.bib literature.bib
	rm -f gitid.tex
	rm -f *~

tar: all
	rm -f $(NAME).tar.gz
	#dirty sed hack ahead to grep files used by latexmk to produce the pdf
	tar cvzhf $(NAME).tar.gz $(NAME).pdf \
  	   `sed -n 's/^[[:space:]]*"\([^/][^"]*\.\(bib\|tex\|eps\|cls\)\)".*$$/\1/p' $(NAME).fdb_latexmk` 

gitid.tex: dummy
	[ -f gitid.tex ] || touch gitid.tex
	echo '$(GITID)' | cmp -s gitid.tex - || echo '$(GITID)' > gitid.tex

dummy: ;
