## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

## set python package name
pkgpythondir = $(pyexecdir)/adaptive_sampling

## compile python files
pkgpython_PYTHON = __init__.py adaptive_sampling.py bayesian_binning.py config.py interface.py policy.py statistics.py visualization.py

## clean python files
clean-local:
	$(RM) *.pyc
