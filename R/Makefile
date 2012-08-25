# Makefile for the R library
R = R
VERSION  = 1.0-1
PACKAGES = adaptive.sampling prombs
TARGETS  = $(addsuffix _$(VERSION).tar.gz,$(PACKAGES))

all: $(TARGETS)

%_$(VERSION).tar.gz: %/NAMESPACE
	$(R) CMD build $*

%/NAMESPACE: %/R/*.R %/src/*.c
	$(R) --vanilla --quiet < roxygen_$*.R

check: $(TARGETS)
	$(R) CMD check $(TARGETS)

install: $(TARGETS)
	$(R) CMD INSTALL $(TARGETS)

clean:
	$(RM) -r $(addsuffix .Rcheck,$(PACKAGES))
	$(RM) $(TARGETS)
