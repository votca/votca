@SET_MAKE@
ACLOCAL_AMFLAGS = -I config
AUTOMAKE_OPTIONS = foreign

SUBDIRS = src include share scripts 

EXTRA_DIST = CHANGELOG NOTICE LICENSE

CLEANFILES = *~
#DISTCLEANFILES =

if HAVE_TXT2TAGS
CHANGELOG: CHANGELOG.t2t
	$(TXT2TAGS) --target txt --infile $< --outfile $@

wiki: CHANGELOG.wiki
%.wiki: %.t2t
	$(TXT2TAGS) --target wiki --infile $< --outfile $@
else
CHANGELOG: CHANGELOG.t2t
	cp $< $@
	echo "For a nicer $@ install txt2tags !"
endif
