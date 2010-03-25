SHELL=/bin/bash
HGID:=$(shell hg parents -R . --template "{node|short}" | sed 's/.*/\\newcommand{\\hgid}{&}/')

NAME=manual
all: $(NAME).pdf

$(NAME).tex: hgid.tex fig_submake functionality_submake reference_submake usage_submake

%.pdf: %.tex
	./latexmk.pl -pdfdvi $<

%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* clean

qclean:
	./latexmk.pl -C $(NAME).tex

clean: qclean fig_subclean functionality_subclean reference_subclean usage_submake
	rm -f $(NAME).fdb_latexmk
	rm -f hgid.tex
	rm -f *~

hgid.tex: dummy
	[ -f hgid.tex ] || touch hgid.tex
	echo '$(HGID)' | cmp -s hgid.tex - || echo '$(HGID)' > hgid.tex

dummy: ;
