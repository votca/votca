all: programs_submake scripts_submake
	$(MAKE) $(MFLAGS) -C xml

%_submake:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $* clean

clean: programs_subclean scripts_subclean
	$(MAKE) $(MFLAGS) -C xml clean
	rm -f *~
