templatedir = ${datadir}/votca/template

template_DATA = csg_template.cc Makefile.template 

install-data-hook:
	cd $(DESTDIR)$(templatedir) && mv Makefile.template Makefile

CLEANFILES = *~

