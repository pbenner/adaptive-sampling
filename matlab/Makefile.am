## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CFLAGS) $(CFLAGS)

MATLAB_LDFLAGS_EXTRA = ../src/.libs/libadaptive-sampling.a

if ENABLE_MINI_GSL
LDFLAGS += $(top_builddir)/libmini-gsl/libmini-gsl.la
endif

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
