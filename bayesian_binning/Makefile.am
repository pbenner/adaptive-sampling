## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## bayesian-binning library
lib_LTLIBRARIES  = libbayesian-binning.la
libbayesian_binning_la_SOURCES = \
	bayesian-binning.c bayesian-binning.h \
	bayesian-binning-test.h bayesian-binning-test.c \
	break-probabilities.c break-probabilities.h \
	datatypes.h \
	effective-counts.c effective-counts.h \
	entropy.c entropy.h \
	interface.c interface.h \
	marginal.c marginal.h \
	model.c model.h \
	model-posterior.c model-posterior.h \
	moment.c moment.h \
	utility.h
libbayesian_binning_la_LIBADD  = $(LIB_PTHREAD)
libbayesian_binning_la_LIBADD += $(top_builddir)/libexception/libexception.la
libbayesian_binning_la_LIBADD += $(top_builddir)/libprombs/libprombs.la
libbayesian_binning_la_LIBADD += $(top_builddir)/libmgs/libmgs.la
libbayesian_binning_la_LDFLAGS = -no-undefined -version-info 0:0:0
# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libbayesian_binning_la_SOURCES = dummy.cxx	

## install headers
nobase_include_HEADERS = bayesian-binning.h

## set python package name
pkgpythondir = $(pyexecdir)/bayesian_binning

## compile python files
pkgpython_PYTHON = __init__.py adaptive_sampling.py bayesian_binning.py config.py interface.py statistics.py visualization.py

## clean python files
clean-local:
	$(RM) *.pyc
