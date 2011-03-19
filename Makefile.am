## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = libexception libprombs bayesian_binning include

dist_bin_SCRIPTS  = bayesian-binning adaptive-sampling
CLEANFILES        = $(dist_bin_SCRIPTS)

man1_MANS         = adaptive-sampling.1 bayesian-binning.1

EXTRA_DIST        = adaptive-sampling.in bayesian-binning.in    \
                    debian/control                              \
                    debian/copyright                            \
                    debian/compat                               \
                    debian/changelog                            \
                    debian/libbayesian-binning0-dev.install     \
                    debian/rules                                \
                    debian/bayesian-binning.install             \
                    debian/libbayesian-binning0.install         \
                    debian/docs                                 \
                    debian/adaptive-sampling.install            \
                    debian/dirs


DB2MAN            = /usr/share/sgml/docbook/stylesheet/xsl/nwalsh/manpages/docbook.xsl
XP                = xsltproc -''-nonet -''-param man.charmap.use.subset "0"

do_subst = sed	-e 's,[@]pythondir[@],$(pythondir),g'

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
