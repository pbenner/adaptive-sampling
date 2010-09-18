[1]dnl -- configure.ac --
dnl
dnl Copyright (C) 2010 Philipp Benner
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

AC_INIT([Bayesian-Toolbox], [1.0], [philipp.benner@mis.mpg.de])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR(bayesian-binning/bayesian-binning.c)
AM_INIT_AUTOMAKE([gnu dist-bzip2])
AC_PROG_CC
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
AC_ISC_POSIX
AC_PROG_INSTALL
AM_PROG_CC_C_O

# ------------------] GENERAL PACKAGE INFO [

AC_SUBST([PACKAGE_SUMMARY], ["template project"])
AC_SUBST([PACKAGE_URL], ["http://www.homepage.ch"])
AC_SUBST([PACKAGE_LICENSE], ["GPL"])
AC_SUBST([PACKAGE_DESCRIPTION], ["This is a template project using autoconf, automake and autoheader."])

# ------------------] MACRO DEFINITIONS [

AC_DEFUN([AC_MY_ARG_WITH],
	[# create help entry and conditionals
		AC_ARG_WITH([$1],
			AC_HELP_STRING([--with-$1],
				[# help string
					build $1 (default is $2)]),
			ac_cv_use_$1=$withval,
			ac_cv_use_$1=$2)
		AC_CACHE_CHECK(whether to use $1, ac_cv_use_$1, ac_cv_use_$1=$2)])

dnl AC_NEW_SUBDIR(subdir_name, default_config)
AC_DEFUN([AC_NEW_SUBDIR],
        [AC_MY_ARG_WITH([$1], [$3])
	AM_CONDITIONAL([COND_$2], [test "$ac_cv_use_$1" = yes])])

# ------------------] LIBRARIES [


dnl HP-UX provides a library calles libhplx includes linux apis.
AC_ARG_VAR(HPLX_LIBRARY_PATH, [Path to the hplx library directory (only for hp-ux systems).])
AC_ARG_VAR(HPLX_INCLUDE_PATH, [Path to the hplx include directory (only for hp-ux systems).])
AC_CHECK_LIB([hplx], [hstrerror],
	[# Add additional include an library directory.
	LDFLAGS="$LDFLAGS -L${HPLX_LIBRARY_PATH:-/usr/local/hplx/lib}"
	CFLAGS="$CFLAGS -I${HPLX_INCLUDE_PATH:-/usr/local/hplx/include}"],,
	[# Additional libraries:
	-lsec])

AC_CHECK_LIB([ucb], [sprintf])
AC_CHECK_LIB([socket], [connect])
dnl On Solaris systems gethostbyname() is in libresolv.
AC_CHECK_LIB([resolv], [gethostbyname],,
	AC_CHECK_LIB([resolv], [hstrerror]))
AC_CHECK_LIB([nsl], [gethostbyname])
AC_CHECK_LIB([c], [dlopen])
dnl Linux distributes dlopen() in libdl.
AC_CHECK_LIB([dl], [dlopen])
dnl Tru64 distributes vasprintf() in -liberty
AC_CHECK_LIB([iberty], [vasprintf])

# ------------------] GETOPT_LONG [

AC_DEFUN([FUNC_GETOPT_LONG],
	[AC_PREREQ(2.49)dnl
		# clean out junk possibly left behind by a previous configuration
		rm -f include/getopt.h
		# Check for getopt_long support
		AC_CHECK_HEADERS([getopt.h])
		AC_CHECK_FUNCS([getopt_long],,
			[# FreeBSD has a gnugetopt library for this
				AC_CHECK_LIB([gnugetopt],[getopt_long],[AC_DEFINE([HAVE_GETOPT_LONG])],
					[# use the GNU replacement
						AC_LIBOBJ(getopt)
						AC_LIBOBJ(getopt1)
						AC_CONFIG_LINKS([include/getopt.h:include/gnugetopt.h])])])])

FUNC_GETOPT_LONG

# ------------------] DEFINES [
#

AH_VERBATIM([_GNU_SOURCE],
[/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE
#endif])

AC_DEFINE_UNQUOTED(TARGET, "$target",
	[System type this program will compile code for (only needed for compilers).])
AC_DEFINE_UNQUOTED(TARGET_CPU, "$target_cpu",
	[System cpu this program will compile code for (only needed for compilers).])
AC_DEFINE_UNQUOTED(TARGET_OS, "$target_os",
	[Operatins system this program will compile code for (only needed for compilers).])
AC_DEFINE_UNQUOTED(TARGET_VENDOR, "$target_vendor",
	[System vendor type this program will compile code for (only needed for compilers).])

AC_DEFINE_UNQUOTED(BUILD, "$build",
	[System type this program was built on.])
AC_DEFINE_UNQUOTED(BUILD_CPU, "$build_cpu",
	[System cpu this program was built on.])
AC_DEFINE_UNQUOTED(BUILD_OS, "$build_os",
	[Operatins system this program was built on.])
AC_DEFINE_UNQUOTED(BUILD_VENDOR, "$build_vendor",
	[System vendor type this program was built on.])

AC_DEFINE_UNQUOTED(HOST, "$host",
	[System type this program will run on.])
AC_DEFINE_UNQUOTED(HOST_CPU, "$host_cpu",
	[System cpu this program will run on.])
AC_DEFINE_UNQUOTED(HOST_OS, "$host_os",
	[Operatins system this program will run on.])
AC_DEFINE_UNQUOTED(HOST_VENDOR, "$host_vendor",
	[System vendor type this program will run on.])

# ------------------] HEADERS [
#

HEADERS_TO_CHECK="stdarg.h dlfcn.h dirent.h netinet/in.h sys/ioctl.h
									arpa/inet.h netdb.h syslog.h signal.h"

AC_CHECK_HEADERS($HEADERS_TO_CHECK,,
	AC_MSG_ERROR([this header is required]))

AC_CHECK_HEADER(net/if.h,,
	AC_MSG_ERROR([this header is required]),
	[[#include <sys/types.h>
		#include <sys/socket.h>]] )

AC_CHECK_HEADER([error.h],
	AC_DEFINE(HAVE_ERROR_H, [], [Needed on OpenBSD systems.]))

# ------------------] FUNCTIONS [
#

AC_CHECK_FUNCS(dprintf)
AC_CHECK_FUNCS(vsyslog)

# ------------------] SUBDIR CONFIGURATION [
#

#AC_NEW_SUBDIR([doc], [DOC], ["no"])

# ------------------] PERL [
#

AC_PATH_PROG(PERL, perl, "no")
AM_CONDITIONAL([COND_PERL], [test "$PERL" != no])

# ------------------] AC_OUTPUT [
#

AC_CONFIG_HEADER([config.h:config.h.in])
AC_CONFIG_FILES([
	Makefile
	version.h
	libexception/Makefile
	bayesian-binning/Makefile
 	])

dnl AC_CONFIG_FILES([updatedd-wrapper/create-config],
dnl 	[chmod 755 updatedd-wrapper/create-config])
dnl AC_CONFIG_FILES([updatedd-wrapper/updatedd-wrapper],
dnl 	[chmod 755 updatedd-wrapper/updatedd-wrapper])
dnl AC_CONFIG_FILES([scripts/ipserv.pl],
dnl 	[chmod 755 scripts/ipserv.pl])
dnl AC_CONFIG_FILES([scripts/interface.pl],
dnl 	[chmod 755 scripts/interface.pl])

AC_OUTPUT

# Local Variables:
# tab-width: 2
# End:
