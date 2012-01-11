## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## bayesian-binning library
lib_LTLIBRARIES  = libadaptive-sampling.la
libadaptive_sampling_la_SOURCES = \
	bayesian-binning.c \
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
	threading.c threading.h \
	utility.h
libadaptive_sampling_la_LIBADD  = $(LIB_PTHREAD)
libadaptive_sampling_la_LIBADD += $(top_builddir)/libexception/libexception.la
libadaptive_sampling_la_LIBADD += $(top_builddir)/libprombs/libprombs.la
libadaptive_sampling_la_LIBADD += $(top_builddir)/libmgs/libmgs.la
libadaptive_sampling_la_LDFLAGS = -no-undefined -version-info 0:0:0
