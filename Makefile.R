# Makefile for the R library
R = R
LIB = adaptive.sampling
MAN = $(LIB)/man
SRCFILES = $(shell ls $(LIB)/R)
SRC = $(SRCFILES:%=$(LIB)/R/%)
VERSION = 1.0-1
TARGET = $(LIB)_$(VERSION).tar.gz

all: $(TARGET)

$(TARGET): $(LIB) $(MAN)
	if [ -d "$$HOME/.usr" ] ; then \
		export PATH="$$HOME/.usr/bin:$$HOME/.cabal/bin/:$$PATH"; \
		export CPATH=$$HOME/.usr/include; \
		export LD_LIBRARY_PATH=$$HOME/.usr/lib:$$LD_LIBRARY_PATH; \
		export LIBRARY_PATH=$$LD_LIBRARY_PATH; \
		export PYTHONPATH=$$HOME/.usr/lib/python2.6/site-packages/; \
		export PERL5LIB=$$HOME/.usr/share/perl; \
		export MANPATH=$$HOME/.usr/share/man:$$MANPATH; \
	fi ; 
	$(R) CMD build $(LIB)

test: $(TARGET)
	$(R) CMD check $^

install: $(TARGET)
	$(R) CMD INSTALL $^

$(MAN): $(SRC)
	$(R) CMD roxygen -d $(LIB)


