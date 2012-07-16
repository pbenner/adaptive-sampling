SHELL		= /bin/sh -e

LATEX		= pdflatex -halt-on-error -shell-escape
TEXI2DVI	= texi2dvi -p
BIBTEX		= bibtex

TEXINPUTS	= .:../share:../../share:

TEX_OUTPUT	= nips2011.pdf

## default build rule
all: $(TEX_OUTPUT) 
	@echo "Usage: make <pdf-output>"
	@echo "       make html"

## latex&co build rules
$(TEX_OUTPUT): $(wildcard *.texinput)
$(TEX_OUTPUT): $(wildcard *.bib)
$(TEX_OUTPUT): $(wildcard example/*.pdf)
$(TEX_OUTPUT): nips2011_macros.sty
$(TEX_OUTPUT): %.pdf: %.tex
# texi2dvi runs latex as often as needed and generates a
# bibliography by invoking bibtex
	TEXINPUTS=$(TEXINPUTS) LATEX="$(LATEX)" BIBTEX=$(BIBTEX) $(TEXI2DVI) $<

## clean rules
clean:
	$(RM) -r $(TEX_OUTPUT:.pdf=)
	$(RM) -r input/auto
	$(RM) -r auto
	$(RM) tmp.inputs
	$(RM) *.table
	$(RM) *.gnuplot
	$(RM) *.aux
	$(RM) *.blg
	$(RM) *.bbl
	$(RM) *.dvi
	$(RM) *.log
	$(RM) *.ps
	$(RM) *.toc
	$(RM) *.out
	$(RM) *.nav
	$(RM) *.snm
	$(RM) *.nlo
	$(RM) *.nls
	$(RM) *.ilg
	$(RM) *.idx
	$(RM) *.ind
	$(RM) *.thm
	$(RM) -f $(ORIGINAL)

distclean: clean
	$(RM) $(TEX_OUTPUT)

.PHONY: all clean distclean html
