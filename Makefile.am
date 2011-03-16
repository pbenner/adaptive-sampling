## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = libexception libprombs bayesian_binning include

dist_bin_SCRIPTS  = bayesian-binning adaptive-sampling
CLEANFILES        = $(dist_bin_SCRIPTS)
EXTRA_DIST        = adaptive-sampling.in bayesian-binning.in

do_subst = sed	-e 's,[@]pythondir[@],$(pythondir),g'

adaptive-sampling: adaptive-sampling.in Makefile
	$(do_subst) < $(srcdir)/adaptive-sampling.in > adaptive-sampling
	chmod +x adaptive-sampling

bayesian-binning: bayesian-binning.in Makefile
	$(do_subst) < $(srcdir)/bayesian-binning.in > bayesian-binning
	chmod +x bayesian-binning
