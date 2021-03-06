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

AC_INIT([Adaptive-Sampling], [1.1], [philipp.benner@mis.mpg.de])
AC_CONFIG_AUX_DIR([auxtools])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([src/main.c])
dnl we have a README.md instead of a README
AM_INIT_AUTOMAKE([gnu dist-bzip2 foreign])
dnl turn of AC_PROG_CC's default behaviour
: ${CFLAGS=""}
AC_PROG_CC
AC_C_CONST
AC_ISC_POSIX
AC_PROG_INSTALL
AM_PROG_CC_C_O
AC_PROG_LN_S
dnl libtool
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL
LT_INIT

dnl ,---------------------------- 
dnl | GENERAL PACKAGE INFO
dnl `----------------------------

AC_SUBST([PACKAGE_SUMMARY], ["Adaptive Sampling"])
AC_SUBST([PACKAGE_URL], ["http://www.mis.mpg.de"])
AC_SUBST([PACKAGE_LICENSE], ["GPL"])
AC_SUBST([PACKAGE_DESCRIPTION], ["Package for adaptive sampling of psychometric functions."])

dnl ,---------------------------- 
dnl | MACRO DEFINITIONS
dnl `----------------------------

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

dnl ,---------------------------- 
dnl | LIBRARIES
dnl `----------------------------

dnl HP-UX provides a library called libhplx which provides linux apis.
AC_ARG_VAR(HPLX_LIBRARY_PATH, [Path to the hplx library directory (only for hp-ux systems).])
AC_ARG_VAR(HPLX_INCLUDE_PATH, [Path to the hplx include directory (only for hp-ux systems).])
AC_CHECK_LIB([hplx], [hstrerror],
	[# Add additional include an library directory.
	LDFLAGS="$LDFLAGS -L${HPLX_LIBRARY_PATH:-/usr/local/hplx/lib}"
	CFLAGS="$CFLAGS -I${HPLX_INCLUDE_PATH:-/usr/local/hplx/include}"],,
	[# Additional libraries:
	-lsec])

AC_CHECK_LIB([ucb], [sprintf])
dnl Tru64 distributes vasprintf() in -liberty
AC_CHECK_LIB([iberty], [vasprintf])
AC_CHECK_LIB([m],[cos],,
  AC_MSG_ERROR([this library is required]))

AC_ARG_ENABLE([pthreads],
	[AS_HELP_STRING([--enable-pthreads],
								[enable pthreads if available @<:@default=yes@:>@])],
	[],
	[enable_pthreads=yes])
AS_IF([test "x$enable_pthreads" = "xyes"],
	[AC_CHECK_LIB([pthread],
								[pthread_create],
								[LIB_PTHREAD="-lpthread"
								 have_pthreads="yes"
								 AC_DEFINE([HAVE_LIB_PTHREAD], [1], [Use pthread library.])],
								 AC_MSG_ERROR([this library is required]))],
	[have_pthreads="no"])

AC_SUBST(LIB_PTHREAD)

dnl ,---------------------------- 
dnl | MINI-GSL LIBRARY
dnl `----------------------------

AC_ARG_ENABLE([mini-gsl],
  [  --enable-mini-gsl    Use mini-gsl instead to dynamically linking to the full gsl. ],
  [case "${enableval}" in
    yes) enable_mini_gsl=yes ;;
    no)  enable_mini_gsl=no ;;
    *) AC_MSG_ERROR([bad value ${enableval} for --enable-mini-gsl]) ;;
  esac],[enable_mini_gsl=no])

AS_IF([test "x$enable_mini_gsl" = "xyes"],[],
	[AC_CHECK_LIB([gslcblas],[cblas_dgemm],,
			AC_MSG_ERROR([this library is required]))
	 AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,
	 		AC_MSG_ERROR([this library is required]))])

AM_CONDITIONAL([ENABLE_MINI_GSL], [test x$enable_mini_gsl = xyes])

dnl ,---------------------------- 
dnl | GETOPT_LONG
dnl `----------------------------

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

dnl ,---------------------------- 
dnl | DEFINES
dnl `----------------------------

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

dnl ,---------------------------- 
dnl | HEADERS
dnl `----------------------------

HEADERS_TO_CHECK="stdarg.h"

AC_CHECK_HEADERS($HEADERS_TO_CHECK,,
	AC_MSG_ERROR([this header is required]))

AC_CHECK_HEADER([errno.h],
	AC_DEFINE(HAVE_ERRNO_H, [], [Needed on OpenBSD systems.]))
AC_CHECK_HEADER([netdb.h],
	AC_DEFINE(HAVE_NETDB_H, [], [Not present on mingw.]))
AC_CHECK_HEADER([syslog.h],
	AC_DEFINE(HAVE_SYSLOG_H, [], [Not present on mingw.]))

dnl ,---------------------------- 
dnl | FUNCTIONS
dnl `----------------------------

AC_CHECK_FUNCS(hstrerror)
AC_CHECK_FUNCS(vsyslog)

dnl ,---------------------------- 
dnl | PRECISION
dnl `----------------------------

AC_ARG_ENABLE([longdouble],
	[AS_HELP_STRING([--enable-longdouble],
								[use long double for probabilities @<:@default=yes@:>@])],
	[],
	[enable_longdouble=yes])
AS_IF([test "x$enable_longdouble" = "xyes"],
	[AC_CHECK_FUNCS([expl logl log1pl], [], [enable_longdouble=no])], [])
AS_IF([test "x$enable_longdouble" = "xyes"],
	[AC_DEFINE([PROB_T], [long double], [Use long double to store log probabilities.])
	 AC_DEFINE([EXP], [expl], [Replace expl with exp])
	 AC_DEFINE([LOG], [logl], [Replace logl with log])
	 AC_DEFINE([LOG1P], [log1pl], [Replace log1pl with log1p])],
	[AC_DEFINE([PROB_T], [double], [Use double to store log probabilities.])
	 AC_DEFINE([EXP], [exp], [Replace expl with exp])
	 AC_DEFINE([LOG], [log], [Replace logl with log])
	 AC_DEFINE([LOG1P], [log1p], [Replace log1pl with log1p])
	 ])

dnl ,---------------------------- 
dnl | GCC_OPTIONS
dnl `----------------------------

# do not reset CFLAGS so we can accept command line options
# CFLAGS=""
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
								[turn on debugging @<:@default=no@:>@])],
	[],
	[enable_debug=no])
AS_IF([test "x$enable_debug" = "xyes"],
  [AX_CFLAGS_GCC_OPTION([-O0])
   AX_CFLAGS_GCC_OPTION([-g])
   AX_CFLAGS_GCC_OPTION([-ggdb])
   AX_CFLAGS_GCC_OPTION([-pg])
	 CFLAGS="$CFLAGS -DDEBUG"],
  [AX_CFLAGS_GCC_OPTION([-O3])
   AX_CFLAGS_GCC_OPTION([-ffast-math])])

AX_CFLAGS_GCC_OPTION([-Wall])
AX_CFLAGS_GCC_OPTION([-Wwrite-strings])
AX_CFLAGS_GCC_OPTION([-Winline])
AX_CFLAGS_GCC_OPTION([-Wstrict-prototypes])
AX_CFLAGS_GCC_OPTION([-Wno-trigraphs])
AX_CFLAGS_GCC_OPTION([-Wtrampolines])
dnl AX_CFLAGS_GCC_OPTION([-ansi])
AX_CFLAGS_GCC_OPTION([-fPIC])
AX_CFLAGS_GCC_OPTION([-fno-nested-functions])

dnl ,---------------------------- 
dnl | MATLAB MEX SUPPORT
dnl `----------------------------

MEX_INIT([MEX],[1],[1],[2])
CHECK_MATLAB

dnl ,---------------------------- 
dnl | SUBDIR CONFIGURATION
dnl `----------------------------

#AC_NEW_SUBDIR([doc], [DOC], ["no"])

dnl ,---------------------------- 
dnl | PERL
dnl `----------------------------

dnl AC_PATH_PROG(PERL, perl, "no")
dnl AM_CONDITIONAL([COND_PERL], [test "$PERL" != no])

dnl ,---------------------------- 
dnl | PYTHON
dnl `----------------------------

AM_PATH_PYTHON(2.5)

dnl ,---------------------------- 
dnl | CYGWIN
dnl `----------------------------

dnl disable library dependency check on cygwin
dnl because of problems with -liberty
AS_IF([test "x$host_os" = "xcygwin"],
  [deplibs_check_method="pass_all"], [])

dnl ,---------------------------- 
dnl | AC_OUTPUT
dnl `----------------------------

AC_CONFIG_HEADER([config.h:config.h.in])
AC_CONFIG_HEADER([probtype.h:probtype.h.in])
AC_CONFIG_LINKS([include/adaptive-sampling/probtype.h:probtype.h])
AC_CONFIG_FILES([
	Makefile
	version.h
	auxtools/Makefile
	data/Makefile
	include/Makefile
	libprombs/Makefile
	libmgs/Makefile
	libexception/Makefile
	libmini-gsl/Makefile
	adaptive_sampling/Makefile
	matlab/Makefile
	math/Makefile
	src/Makefile
 	])

AC_OUTPUT

dnl ,---------------------------- 
dnl | Print a summary
dnl `----------------------------

AC_MSG_NOTICE_BOX([Configure summary

System type............... ${host}
Build with debug symbols.. ${enable_debug}
Use long double........... ${enable_longdouble}
Matlab support............ ${have_matlab}
PThread support........... ${have_pthreads}
Use mini-gsl.............. ${enable_mini_gsl}
])

dnl Local Variables:
dnl tab-width: 2
dnl End:
