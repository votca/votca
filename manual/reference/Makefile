SHELL=/bin/bash
all: programs_submake xml_submake calculators_submake

%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* clean

dummy: ;

clean: programs_subclean xml_subclean calculators_subclean
	rm -f *~ .log *.backup *.aux
