## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

MATLAB_LINK_EXTRA = -Wl,--no-undefined -Wl,-dn,../src/.libs/libadaptive-sampling.a -static-libgcc -Wl,-dn,-lgsl,-dy,-lmwblas

all: example.$(MEXEXT)

example.o: example.c
	$(CC) $(DEFS) $(CFLAGS) $(INCLUDES) -I$(top_srcdir) $(MATLAB_FLAGS) -c -o $@ $<
example.$(MEXEXT): example.o
	$(CC) -o $@ $< $(MATLAB_LIB) $(MATLAB_LINK) $(MATLAB_LINK_EXTRA) $(LIB_PTHREAD)

clean-local:
	$(RM) example.o
	$(RM) example.$(MEXEXT)
