.PHONY: clean

%_submake:
	$(MAKE) $(MFLAGS) -C $*

%_subclean:
	$(MAKE) $(MFLAGS) -C $* clean

all: fig_submake

clean: fig_subclean
	rm -f *~
