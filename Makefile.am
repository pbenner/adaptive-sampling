## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = include libexception libprombs libmgs
if ENABLE_MINI_GSL
SUBDIRS += libmini-gsl
endif
SUBDIRS += adaptive_sampling src

dist_bin_SCRIPTS  = bayesian-binning adaptive-sampling
CLEANFILES        = $(dist_bin_SCRIPTS)

man1_MANS         = adaptive-sampling.1 bayesian-binning.1

EXTRA_DIST        = adaptive-sampling.1 bayesian-binning.1      \
                    adaptive-sampling.in bayesian-binning.in    \
                    adaptive-sampling.xml bayesian-binning.xml  \
                    data/binning/*.cfg                          \
                    data/binning/*.py                           \
                    data/binning/*.pdf                          \
                    data/binning-spikes/*.cfg                   \
                    data/binning-spikes/*.py                    \
                    data/binning-spikes/*.pdf                   \
                    data/binning-spikes/*.result*               \
                    data/sampling/*.cfg                         \
                    data/sampling/*.py                          \
                    data/sampling/*.pdf                         \
                    data/sampling/*.result*                     \
                    debian/control                              \
                    debian/copyright                            \
                    debian/compat                               \
                    debian/changelog                            \
                    debian/libbayesian-binning0-dev.install     \
                    debian/rules                                \
                    debian/bayesian-binning.examples            \
                    debian/bayesian-binning.install             \
                    debian/bayesian-binning.README.Debian       \
                    debian/libbayesian-binning0.install         \
                    debian/docs                                 \
                    debian/adaptive-sampling.examples           \
                    debian/adaptive-sampling.install            \
                    debian/adaptive-sampling.README.Debian

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

