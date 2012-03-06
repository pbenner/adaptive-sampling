## Process this file with automake to produce Makefile.in

## bayesian-binning library
lib_LTLIBRARIES  = libadaptive-sampling.la
libadaptive_sampling_la_SOURCES = \
	bayesian-binning.c \
	bayesian-binning-test.h bayesian-binning-test.c \
	break-probabilities.c break-probabilities.h \
	datatypes.h \
	effective-counts.c effective-counts.h \
	interface.c interface.h \
	marginal.c marginal.h \
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
