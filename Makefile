SHELL=/bin/bash
#the next line is used by the buildutil !
VER=1.2.1
HGID:=$(shell hg parents -R . --template "{node|short}" | sed 's/.*/\\newcommand{\\hgid}{${VER} (&)}/')
LATEXMK=scripts/latexmk.pl
LATEXMKOPTS=-e '$$latex=q/latex --halt-on-error %O %S/'

NAME=manual
all: $(NAME).pdf
dvi: $(NAME).dvi
ps: $(NAME).ps

$(NAME).tex: hgid.tex fig_submake functionality_submake reference_submake usage_submake

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
	echo "$(NAME).pdf is ready to use, no need to install it."

clean: qclean fig_subclean functionality_subclean reference_subclean usage_subclean
	rm -f $(NAME).fdb_latexmk $(NAME).brf
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

upload-pdf: manual.pdf
	googlesites_upload.py -d "/Documentation" -a manual.pdf

dummy: ;
