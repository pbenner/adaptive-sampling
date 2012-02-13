dnl Copyright (C) 2012 Philipp Benner
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

AC_DEFUN([CHECK_MATLAB],
[
AC_ARG_WITH(matlab,
	[  --with-matlab=DIR	the directory where Matlab is installed ],
	[if test -n "${withval}" -a ! "${withval}" = "yes"; then
		MATLAB=${withval};
	 else
		AC_PATH_PROG([MATLAB_PROG], [matlab])
		if test -n "${MATLAB_PROG}"; then
		   	MATLAB=$($MATLAB_PROG -e | grep '^MATLAB=' | sed s,MATLAB=,,);
		else
			AC_MSG_ERROR([You must specify the matlab directory with --with-matlab=PATH])
		fi;
	 fi],
	MATLAB=)
AC_MSG_CHECKING(for Matlab)

dnl ,---------------------------- 
dnl | Source matlab variables
dnl `----------------------------
if test -n "${MATLAB}"
then
	eval $(./auxtools/matlab.sh ${MATLAB})
	AC_SUBST(MATLAB_RPATH)
	AC_SUBST(MATLAB_CFLAGS)
	AC_SUBST(MATLAB_LIBS)
	AC_SUBST(MATLAB_LDFLAGS)
	AC_SUBST(MATLAB_LDEXTENSION)
	AM_CONDITIONAL(HAVE_MATLAB, test "xyes" = "xyes" )
	have_matlab=yes
else
	AM_CONDITIONAL(HAVE_MATLAB, test "xno" = "xyes" )
	have_matlab=no
fi
AC_MSG_RESULT($have_matlab)
])
