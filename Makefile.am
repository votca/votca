ACLOCAL_AMFLAGS = -I config

#script should be before src because help2t2t to build manpages
SUBDIRS = scripts src include share

EXTRA_DIST = CHANGELOG.t2t CHANGELOG NOTICE LICENSE

CLEANFILES = *~
DISTCLEANFILES =  CHANGELOG

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
