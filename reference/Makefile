SHELL=/bin/bash
all: check csgid.tex programs_submake scripts_submake xml_submake 

upload: scripts_upload xml_upload programs_upload


check:
	@if ! which csg_call &>/dev/null; then\
	  echo '#################################';\
	  echo '#                               #';\
	  echo '#  Could not find csg_call      #';\
	  echo '#  Please source VOTCARC first  #';\
	  echo '#                               #';\
	  echo '#################################';\
	  exit 1;\
	fi

%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_upload:
	$(MAKE) $(MFLAGS) -C $* upload


%_subclean:
	$(MAKE) $(MFLAGS) -C $* clean

CSGID:=$(shell csg_call --help | sed -n 's/.*\(hg\|git\)id: \(.*\)/\\newcommand{\\csgid}{\2}/p')
csgid.tex: check dummy
	[ -f csgid.tex ] || touch csgid.tex
	echo '$(CSGID)' | cmp -s csgid.tex - || echo '$(CSGID)' > csgid.tex

dummy: ;

clean: programs_subclean scripts_subclean xml_subclean
	rm -f *~ csgid.tex
