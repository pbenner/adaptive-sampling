## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

EXTRA_DIST = \
	README \
	interface.h \
	binningPosterior.c \
	samplingUtility.c \
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

all: binningPosterior$(MATLAB_LDEXTENSION) samplingUtility$(MATLAB_LDEXTENSION)

binningPosterior.o: binningPosterior.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
binningPosterior$(MATLAB_LDEXTENSION): binningPosterior.o
	$(CC) -o $@ $< $(LDFLAGS) $(MATLAB_RPATH) $(MATLAB_LDFLAGS) $(MATLAB_LDFLAGS_EXTRA) $(MATLAB_LIBS) $(LIBS) $(LIB_PTHREAD)
samplingUtility.o: samplingUtility.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
samplingUtility$(MATLAB_LDEXTENSION): samplingUtility.o
	$(CC) -o $@ $< $(LDFLAGS) $(MATLAB_RPATH) $(MATLAB_LDFLAGS) $(MATLAB_LDFLAGS_EXTRA) $(MATLAB_LIBS) $(LIBS) $(LIB_PTHREAD)

tar:
	mkdir adaptive-sampling-matlab
	cp  *.m *.mex* adaptive-sampling-matlab
	tar -cjvf adaptive-sampling-matlab.tar.bz2 adaptive-sampling-matlab
	$(RM) -r adaptive-sampling-matlab

clean-local:
	$(RM) binningPosterior.o
	$(RM) binningPosterior$(MATLAB_LDEXTENSION)
	$(RM) samplingUtility.o
	$(RM) samplingUtility$(MATLAB_LDEXTENSION)
