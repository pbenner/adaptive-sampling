# Makefile for the R library
R          = R
RFLAGS     = --no-multiarch
VERSION    = 1.0-1
PACKAGES   = adaptive.sampling prombs
TARGETS    = $(addsuffix _$(VERSION).tar.gz,$(PACKAGES))
NAMESPACES = $(addsuffix /NAMESPACE,$(PACKAGES))

all: $(TARGETS)

$(TARGETS): %_$(VERSION).tar.gz: %/NAMESPACE
	$(R) CMD build $(RFLAGS) $*

$(NAMESPACES): %/NAMESPACE: $(wildcard %/R/*.R) $(wildcard %/src/*.c)
	$(R) --vanilla --quiet < roxygen_$*.R

check: $(TARGETS)
	$(R) CMD check $(RFLAGS) $(TARGETS)

install: $(TARGETS)
	$(R) CMD INSTALL $(RFLAGS) $(TARGETS)

clean:
	$(RM) -r $(addsuffix .Rcheck,$(PACKAGES))
	$(RM) $(TARGETS)
	$(RM) $(NAMESPACES)
