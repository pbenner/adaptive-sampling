## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = libexception bayesian-binning

## relocate distdir
#distdir = ../$(PACKAGE)-$(VERSION)

## include debian files
EXTRA_DIST = debian/compat debian/control debian/copyright debian/postinst debian/rules \
	     debian/updatedd.files debian/updatedd-doc.files debian/libupdatedd0.files \
	     debian/libupdatedd0-dev.files debian/changelog
