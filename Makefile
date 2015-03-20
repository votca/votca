SHELL=/bin/bash
#the next line is used by the buildutil !
VER=1.3-dev
GITID_PLAIN:=$(shell git rev-parse --short HEAD 2> /dev/null || hg parents -R . --template "{node|short}")
GITID:=$(shell echo $(GITID_PLAIN) | sed 's/.*/\\newcommand{\\gitid}{${VER} (&)}/')
LATEXMK=scripts/latexmk.pl
LATEXMKOPTS=-e '$$latex=q/latex --halt-on-error %O %S/'

NAME=manual
all: $(NAME).pdf
dvi: $(NAME).dvi
ps: $(NAME).ps

$(NAME).tex: gitid.tex fig_submake functionality_submake reference_submake usage_submake

%.dvi: %.tex dummy
	@#rm target if latexmk failed, worked better than DELETE_ON_ERROR
	$(LATEXMK) $(LATEXMKOPTS) -dvi $< || rm -f $@
	@#rm has exit code 0
	@[ -f $@ ]

%.pdf: %.dvi
	dvipdf $< $*_shadow.pdf
	mv $*_shadow.pdf $@

%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* clean

qclean:
	$(LATEXMK) -C $(NAME).tex

install: all
	@echo "$(NAME).pdf is ready to use, no need to install it."

clean: qclean fig_subclean functionality_subclean reference_subclean usage_subclean
	rm -f $(NAME).fdb_latexmk $(NAME).brf
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

upload-pdf: manual.pdf
	googlesites_upload.py -d "/Documentation" -a manual.pdf

dummy: ;
