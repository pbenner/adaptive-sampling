## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## bayesian-binning library
lib_LTLIBRARIES  = libbayesian-binning.la
libbayesian_binning_la_SOURCES = bayesian-binning.c bayesian-binning.h interface.c interface.h
libbayesian_binning_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libbayesian_binning_la_LIBADD += $(top_builddir)/libprombs/libprombs.la

## compile python files
pkgpython_PYTHON = config.py interface.py statistics.py visualization.py

## create symbolic link
all-local: libbayesian-binning.so

libbayesian-binning.so:
	$(LN_S) .libs/libbayesian-binning.so libbayesian-binning.so

clean-local:
	$(RM) libbayesian-binning.so
	$(RM) __init__.pyc
	$(RM) config.pyc
	$(RM) interface.pyc
	$(RM) visualization.pyc
	$(RM) statistics.pyc
