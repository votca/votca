all: programs_submake scripts_submake xml_submake

%_submake:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $* clean

clean: programs_subclean scripts_subclean xml_subclean
	rm -f *~
