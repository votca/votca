SHELL=/bin/bash
all: check hgid.tex programs_submake scripts_submake xml_submake 

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

HGID:=$(shell csg_call --help | sed -n 's/.*hgid: \(.*\)/\\newcommand{\\refhgid}{\1}/p')
hgid.tex: check dummy
	[ -f hgid.tex ] || touch hgid.tex
	echo '$(HGID)' | cmp -s hgid.tex - || echo '$(HGID)' > hgid.tex

dummy: ;

clean: programs_subclean scripts_subclean xml_subclean
	rm -f *~ hgid.tex
