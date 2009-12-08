.PHONY: clean

%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* clean

all: abb_submake

clean: abb_subclean
	rm -f *~
