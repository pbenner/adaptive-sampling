## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## bayesian-binning library
noinst_LTLIBRARIES  = libbayesian-binning.la
libbayesian_binning_la_SOURCES = bayesian-binning.c bayesian-binning.h
libbayesian_binning_la_LIBADD  = $(top_builddir)/libexception/libexception.la
libbayesian_binning_la_LIBADD += $(top_builddir)/libprombs/libprombs.la

## main program
bin_PROGRAMS = bayesian-binning
bayesian_binning_SOURCES = main.c
bayesian_binning_LDADD   = $(LIBOBJS)
bayesian_binning_LDADD  += libbayesian-binning.la
bayesian_binning_LDADD  += $(top_builddir)/libexception/libexception.la
bayesian_binning_LDADD  += $(top_builddir)/libprombs/libprombs.la

## swig
BUILT_SOURCES    = interface_wrap.c
SWIG_SOURCES     = interface.i interface.c

pkgpython_PYTHON        = interface.py
pkgpyexec_LTLIBRARIES   = _interface.la
_interface_la_SOURCES   = interface_wrap.c $(SWIG_SOURCES)
_interface_la_CFLAGS    = $(AX_SWIG_PYTHON_CPPFLAGS)
_interface_la_LDFLAGS   = -module
_interface_la_LIBADD     = libbayesian-binning.la

interface_wrap.c: $(SWIG_SOURCES)
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -o $@ $<

all-local: _interface.so

_interface.so:
	$(LN_S) .libs/_interface.so _interface.so

clean-local:
	$(RM) __init__.pyc
	$(RM) config.pyc
	$(RM) visualization.pyc
	$(RM) statistics.pyc
	$(RM) _interface.so
	$(RM) interface.py
	$(RM) interface.pyc
	$(RM) interface_wrap.c

## manuals
man1_MANS = bayesian-binning.1

## include extra files
EXTRA_DIST = bayesian-binning.1
