## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## bayesian-binning library
lib_LTLIBRARIES  = libadaptive-sampling.la
libadaptive_sampling_la_SOURCES = \
	break-probabilities.c break-probabilities.h \
	datatypes.h \
	density.c density.h \
	effective-counts.c effective-counts.h \
	interface.c interface.h \
	main.c \
	main-test.h main-test.c \
	model.c model.h \
	model-posterior.c model-posterior.h \
	moment.c moment.h \
	threading.c threading.h \
	tools.h \
	utility.c utility.h
libadaptive_sampling_la_LIBADD  = $(LIB_PTHREAD)
libadaptive_sampling_la_LIBADD += $(top_builddir)/libexception/libexception.la
libadaptive_sampling_la_LIBADD += $(top_builddir)/libprombs/libprombs.la
libadaptive_sampling_la_LIBADD += $(top_builddir)/libmgs/libmgs.la
if ENABLE_MINI_GSL
libadaptive_sampling_la_LIBADD += $(top_builddir)/libmini-gsl/libmini-gsl.la
endif
libadaptive_sampling_la_LDFLAGS = -no-undefined -version-info 0:0:0
