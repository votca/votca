SUBDIRS =  detail program_options smart_ptr

AM_CPPFLAGS = -I$(srcdir)/..

lib_LTLIBRARIES = libvotca_boost.la

libvotca_boost_la_SOURCES = \
program_options/src/cmdline.cpp program_options/src/config_file.cpp \
program_options/src/convert.cpp \
program_options/src/options_description.cpp \
program_options/src/parsers.cpp program_options/src/positional_options.cpp \
program_options/src/split.cpp program_options/src/utf8_codecvt_facet.cpp \
program_options/src/value_semantic.cpp \
program_options/src/variables_map.cpp program_options/src/winmain.cpp \
smart_ptr/src/sp_collector.cpp smart_ptr/src/sp_debug_hooks.cpp 

