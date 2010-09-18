## Process this file with automake to produce Makefile.in

## common CFLAGS
AM_CFLAGS = -Wall -Wwrite-strings -Wstrict-prototypes -Wno-trigraphs -fPIC

## main program
bin_PROGRAMS = bayesian-binning
bayesian_binning_SOURCES = bayesian-binning.c bayesian-binning.h
bayesian_binning_LDADD = $(LIBOBJS) $(top_builddir)/libexception/libexception.la
bayesian_binning_CFLAGS = -DLIBPATH="\"$(pkglibdir)\""

INCLUDES = -I$(top_srcdir)/include

## manuals
man1_MANS = bayesian-binning.1

## include extra files
EXTRA_DIST = bayesian-binning.1
