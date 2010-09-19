## Process this file with automake to produce Makefile.in

## common CFLAGS
##AM_CFLAGS = -Wall -Wwrite-strings -Wstrict-prototypes -Wno-trigraphs -fPIC

## main program
bin_PROGRAMS = bayesian-binning
bayesian_binning_SOURCES = bayesian-binning.c bayesian-binning.h
bayesian_binning_LDADD = $(LIBOBJS) $(top_builddir)/libexception/libexception.la
bayesian_binning_CFLAGS = -DLIBPATH="\"$(pkglibdir)\""

INCLUDES = -I$(top_srcdir)/include

## swig
BUILT_SOURCES = $(srcdir)/example_wrap.cpp
EXAMPLE_SOURCES  = example.i

pkgpython_PYTHON = example.py
pkgpyexec_LTLIBRARIES = _example.la
_example_la_SOURCES   = $(srcdir)/example_wrap.cpp $(EXAMPLE_SOURCES)
_example_la_CFLAGS    = $(SWIG_PYTHON_CFLAGS)
_example_la_LDFLAGS   = -module

$(srcdir)/example_wrap.cpp: $(EXAMPLE_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src -o $@ $<

## manuals
man1_MANS = bayesian-binning.1

## include extra files
EXTRA_DIST = bayesian-binning.1
