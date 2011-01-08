EXTRA_DIST = Doxyfile

if HAVE_DOXYGEN
html: Doxyfile
	$(DOXYGEN) Doxyfile

clean-local:
	rm -rf html

endif

