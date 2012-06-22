## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

EXTRA_DIST = \
	adaptive_sampling.c \
	adsamp.m \
	categorial_sample.m \
	count_statistic.m \
	default_alpha.m \
	default_gamma.m \
	default_options.m \
	demo.m \
	example1.m \
	example2.m \
	example3.m \
	getNumberOfCores.m \
	plot_binning.m \
	plot_marginal.m

COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(AM_CPPFLAGS) $(CFLAGS) $(AM_CFLAGS)

MATLAB_LDFLAGS_EXTRA = ../src/.libs/libadaptive-sampling.a

all: adaptive_sampling$(MATLAB_LDEXTENSION)

adaptive_sampling.o: adaptive_sampling.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
adaptive_sampling$(MATLAB_LDEXTENSION): adaptive_sampling.o
	$(CC) -o $@ $< $(LDFLAGS) $(MATLAB_RPATH) $(MATLAB_LDFLAGS) $(MATLAB_LDFLAGS_EXTRA) $(MATLAB_LIBS) $(LIBS) $(LIB_PTHREAD)

tar:
	mkdir adaptive-sampling-matlab
	cp  *.m *.mex* adaptive-sampling-matlab
	tar -cjvf adaptive-sampling-matlab.tar.bz2 adaptive-sampling-matlab
	$(RM) -r adaptive-sampling-matlab

clean-local:
	$(RM) adaptive_sampling.o
	$(RM) adaptive_sampling$(MATLAB_LDEXTENSION)
