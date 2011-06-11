## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## bayesian-binning library
lib_LTLIBRARIES  = libbayesian-binning.la
libbayesian_binning_la_SOURCES = bayesian-binning.c bayesian-binning.h interface.c interface.h
libbayesian_binning_la_LIBADD  = $(LIB_PTHREAD)
libbayesian_binning_la_LIBADD += $(top_builddir)/libexception/libexception.la
libbayesian_binning_la_LIBADD += $(top_builddir)/libprombs/libprombs.la
libbayesian_binning_la_LDFLAGS = -no-undefined -version-info 0:0:0

## install headers
nobase_include_HEADERS = bayesian-binning.h

## set python package name
pkgpythondir = $(pyexecdir)/bayesian_binning

## compile python files
pkgpython_PYTHON = __init__.py config.py binning.py interface.py sampling.py statistics.py visualization.py

## clean python files
clean-local:
	$(RM) *.pyc
