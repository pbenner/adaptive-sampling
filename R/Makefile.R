# Makefile for the R library
R = R
VERSION = 1.0-1
TARGET1 = adaptive.sampling
TARGET2 = prombs

all: $(TARGET1)_$(VERSION).tar.gz $(TARGET2)_$(VERSION).tar.gz

$(TARGET1)_$(VERSION).tar.gz: $(TARGET1) $(TARGET1)/NAMESPACE $(TARGET1)/man
	$(R) CMD build $(TARGET1)
$(TARGET2)_$(VERSION).tar.gz: $(TARGET2) $(TARGET2)/NAMESPACE $(TARGET2)/man
	$(R) CMD build $(TARGET2)

$(TARGET1)/NAMESPACE $(TARGET1)/man: $(wildcard $(TARGET1)/R/*.R)
	$(R) --vanilla --quiet < roxygen_adaptive.sampling.R
$(TARGET2)/NAMESPACE $(TARGET2)/man: $(wildcard $(TARGET2)/R/*.R)
	$(R) --vanilla --quiet < roxygen_prombs.R

test: $(TARGET1)_$(VERSION).tar.gz $(TARGET2)_$(VERSION).tar.gz
	$(R) CMD check $(TARGET1)
	$(R) CMD check $(TARGET2)

install: $(TARGET1)_$(VERSION).tar.gz $(TARGET2)_$(VERSION).tar.gz
	$(R) CMD INSTALL $(TARGET1)
	$(R) CMD INSTALL $(TARGET2)
