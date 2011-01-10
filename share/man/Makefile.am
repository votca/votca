EXTRA_DIST = votca-csg.t2t

if HAVE_TXT2TAGS
man7_MANS= votca-csg.man

%.man: %.t2t
	${TXT2TAGS} -t man -i $< -o $@
else
man7_MANS=
endif

CLEANFILES = $(man7_MANS)
