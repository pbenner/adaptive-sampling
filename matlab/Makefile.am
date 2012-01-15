## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)

MATLAB_LINK_EXTRA = -Wl,--no-undefined -Wl,-dn,../src/.libs/libadaptive-sampling.a -static-libgcc -Wl,-dn,-lgsl,-dy,-lmwblas

all: adaptive_sampling.$(MEXEXT)

adaptive_sampling.o: adaptive_sampling.c
	$(CC) $(DEFS) $(CFLAGS) $(INCLUDES) $(MATLAB_FLAGS) -c -o $@ $<
adaptive_sampling.$(MEXEXT): adaptive_sampling.o
	$(CC) -o $@ $< $(MATLAB_LIB) $(MATLAB_LINK) $(MATLAB_LINK_EXTRA) $(LIB_PTHREAD)

clean-local:
	$(RM) adaptive_sampling.o
	$(RM) adaptive_sampling.$(MEXEXT)
