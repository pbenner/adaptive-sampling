## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = include libexception libprombs libmgs
if ENABLE_MINI_GSL
SUBDIRS += libmini-gsl
endif
SUBDIRS += adaptive_sampling src

DIST_SUBDIRS = auxtools data include libexception libprombs libmgs libmini-gsl adaptive_sampling src matlab math

dist_bin_SCRIPTS  = bayesian-binning adaptive-sampling
CLEANFILES        = $(dist_bin_SCRIPTS)

man1_MANS         = adaptive-sampling.1 bayesian-binning.1

EXTRA_DIST        = R/adaptive.sampling/R/*.R                   \
                    R/adaptive.sampling/src/*.c                 \
                    R/adaptive.sampling/demo/*.R                \
                    R/adaptive.sampling/demo/00Index            \
                    R/adaptive.sampling/src/Makevars            \
                    R/adaptive.sampling/DESCRIPTION             \
                    R/multibin.summation/R/*.R                  \
                    R/multibin.summation/src/*.c                \
                    R/multibin.summation/demo/*.R               \
                    R/multibin.summation/demo/00Index           \
                    R/multibin.summation/src/Makevars           \
                    R/multibin.summation/DESCRIPTION            \
                    R/Makefile                                  \
                    R/README                                    \
                    R/roxygen_adaptive.sampling.R               \
                    R/roxygen_multibin.summation.R              \
                    adaptive-sampling.1 bayesian-binning.1      \
                    adaptive-sampling.in bayesian-binning.in    \
                    adaptive-sampling.xml bayesian-binning.xml  \
                    debian/control                              \
                    debian/copyright                            \
                    debian/compat                               \
                    debian/changelog                            \
                    debian/libadaptive-sampling0-dev.install    \
                    debian/libadaptive-sampling0.install        \
                    debian/rules                                \
                    debian/docs                                 \
                    debian/adaptive-sampling.examples           \
                    debian/adaptive-sampling.install            \
                    debian/adaptive-sampling.README.Debian      \
                    debian/source/*

DB2MAN            = /usr/share/sgml/docbook/stylesheet/xsl/nwalsh/manpages/docbook.xsl
XP                = xsltproc -''-nonet -''-param man.charmap.use.subset "0"

do_subst          = sed -e 's,[@]pythondir[@],$(pythondir),g'

adaptive-sampling: adaptive-sampling.in Makefile
	$(do_subst) < $(srcdir)/adaptive-sampling.in > adaptive-sampling
	chmod +x adaptive-sampling

bayesian-binning: bayesian-binning.in Makefile
	$(do_subst) < $(srcdir)/bayesian-binning.in > bayesian-binning
	chmod +x bayesian-binning

adaptive-sampling.1: adaptive-sampling.xml
	if [ -f $(DB2MAN) ]; then   \
		$(XP) $(DB2MAN) $<; \
	fi

bayesian-binning.1: bayesian-binning.xml
	if [ -f $(DB2MAN) ]; then   \
		$(XP) $(DB2MAN) $<; \
	fi

all-local:
	@echo
	@echo "You might want to add \`backend : GTKAgg' to"
	@echo "your .matplotlib/matplotlibrc"
	@echo

