## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

EXTRA_DIST = \
	README \
	interface.h \
	interface.c \
	binningPosterior.c \
	samplingUtility.c \
	prombs.c \
	prombsExtended.c \
	prombsTest.m \
	prombsExtendedTest.m \
	adsamp.m \
	categorial_sample.m \
	count_statistic.m \
	default_alpha.m \
	default_gamma.m \
	default_options.m \
	example0.m \
	example1.m \
	example2.m \
	example3.m \
	getNumberOfCores.m \
	plot_binning.m \
	plot_density.m

COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(AM_CPPFLAGS) $(CFLAGS) $(AM_CFLAGS)

MATLAB_LDFLAGS_SAMPLING = ../src/.libs/libadaptive-sampling.a
MATLAB_LDFLAGS_PROMBS   = ../libprombs/.libs/libprombs.a

all: binningPosterior$(MATLAB_LDEXTENSION) samplingUtility$(MATLAB_LDEXTENSION) prombs$(MATLAB_LDEXTENSION) prombsExtended$(MATLAB_LDEXTENSION)

interface.o: interface.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
binningPosterior.o: binningPosterior.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
binningPosterior$(MATLAB_LDEXTENSION): binningPosterior.o interface.o
	$(CC) -o $@ $+ $(LDFLAGS) $(MATLAB_RPATH) $(MATLAB_LDFLAGS) $(MATLAB_LDFLAGS_SAMPLING) $(MATLAB_LIBS) $(LIBS) $(LIB_PTHREAD)
samplingUtility.o: samplingUtility.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
samplingUtility$(MATLAB_LDEXTENSION): samplingUtility.o interface.o
	$(CC) -o $@ $+ $(LDFLAGS) $(MATLAB_RPATH) $(MATLAB_LDFLAGS) $(MATLAB_LDFLAGS_SAMPLING) $(MATLAB_LIBS) $(LIBS) $(LIB_PTHREAD)
prombs.o: prombs.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
prombs$(MATLAB_LDEXTENSION): prombs.o interface.o
	$(CC) -o $@ $+ $(LDFLAGS) $(MATLAB_RPATH) $(MATLAB_LDFLAGS) $(MATLAB_LDFLAGS_PROMBS) $(MATLAB_LIBS) $(LIBS)
prombsExtended.o: prombsExtended.c
	$(COMPILE) $(MATLAB_CFLAGS) -c -o $@ $<
prombsExtended$(MATLAB_LDEXTENSION): prombsExtended.o interface.o
	$(CC) -o $@ $+ $(LDFLAGS) $(MATLAB_RPATH) $(MATLAB_LDFLAGS) $(MATLAB_LDFLAGS_PROMBS) $(MATLAB_LIBS) $(LIBS)

tar:
	mkdir adaptive-sampling-matlab
	cp  *.m ../matlab-bin/*.mex* adaptive-sampling-matlab
	tar -cjvf adaptive-sampling-matlab.tar.bz2 adaptive-sampling-matlab
	$(RM) -r adaptive-sampling-matlab

clean-local:
	$(RM) interface.o
	$(RM) binningPosterior.o
	$(RM) binningPosterior$(MATLAB_LDEXTENSION)
	$(RM) samplingUtility.o
	$(RM) samplingUtility$(MATLAB_LDEXTENSION)
	$(RM) prombs.o
	$(RM) prombs$(MATLAB_LDEXTENSION)
	$(RM) prombsExtended.o
	$(RM) prombsExtended$(MATLAB_LDEXTENSION)
