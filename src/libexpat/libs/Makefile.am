libvotca_expat_la_CPPFLAGS = -I$(srcdir)/../include -DHAVE_EXPAT_CONFIG_H

libvotca_expat_la_LDFLAGS = -no-undefined

lib_LTLIBRARIES = libvotca_expat.la

libvotca_expat_la_SOURCES = xmlparse.c xmlrole.c xmltok.c

EXTRA_HEADERS = xmltok_impl.c xmltok_ns.c

install-exec-hook:
if NO_LA_FILES
	rm -f $(DESTDIR)$(libdir)/libvotca_expat.la
endif

