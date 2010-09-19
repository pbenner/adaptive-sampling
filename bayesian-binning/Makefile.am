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
BUILT_SOURCES    = example_wrap.c
EXAMPLE_SOURCES  = example.i example.c

pkgpython_PYTHON      = example.py
pkgpyexec_LTLIBRARIES = _example.la
_example_la_SOURCES   = example_wrap.c $(EXAMPLE_SOURCES)
_example_la_CFLAGS    = $(AX_SWIG_PYTHON_CPPFLAGS)
_example_la_LDFLAGS   = -module

example_wrap.c: $(EXAMPLE_SOURCES)
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -o $@ $<

clean-local:
	$(RM) example_wrap.c example.py example.pyc

## manuals
man1_MANS = bayesian-binning.1

## include extra files
EXTRA_DIST = bayesian-binning.1
