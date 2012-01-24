## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
        $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)

all: adaptive_sampling.$(MEXEXT)

adaptive_sampling.o: adaptive_sampling.c
	$(COMPILE) $(DEFS) $(CFLAGS) $(INCLUDES) $(MATLAB_FLAGS) -c -o $@ $<
adaptive_sampling.$(MEXEXT): adaptive_sampling.o
	$(COMPILE) -o $@ $< $(LDFLAGS) $(MATLAB_LIB) $(MATLAB_LINK) $(MATLAB_LINK_EXTRA) $(LIB_PTHREAD)

tar:
	mkdir adaptive-sampling-matlab
	cp  *.m *.mex* adaptive-sampling-matlab
	tar -cjvf adaptive-sampling-matlab.tar.bz2 adaptive-sampling-matlab
	$(RM) -r adaptive-sampling-matlab

clean-local:
	$(RM) adaptive_sampling.o
	$(RM) adaptive_sampling.$(MEXEXT)
