SHELL=/bin/bash
all: check programs_submake scripts_submake xml_submake

check:
	@if ! which csg_property &>/dev/null; then\
	  echo '#################################';\
	  echo '#                               #';\
	  echo '#  Could not find csg_property  #';\
	  echo '#  Please source VOTCARC first  #';\
	  echo '#                               #';\
	  echo '#################################';\
	  exit 1;\
	fi

%_submake:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* -f Makefile.in bootstrap
	$(MAKE) $(MFLAGS) -C $* clean

clean: programs_subclean scripts_subclean xml_subclean
	rm -f *~
