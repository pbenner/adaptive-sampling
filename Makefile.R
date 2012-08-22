# Makefile for the R library
R = R
LIB = adaptive.sampling
MAN = $(LIB)/man
SRCFILES = $(shell ls $(LIB)/R)
SRC = $(SRCFILES:%=$(LIB)/R/%)
VERSION = 1.0-1
TARGET = $(LIB)_$(VERSION).tar.gz

all: $(TARGET)

$(TARGET): $(LIB) $(LIB)/NAMESPACE $(MAN)
	$(R) CMD build $(LIB)

test: $(TARGET)
	$(R) CMD check $^

install: $(TARGET)
	$(R) CMD INSTALL $^

$(LIB)/NAMESPACE $(MAN): $(SRC)
	$(R) --vanilla --quiet < roxygen.R
