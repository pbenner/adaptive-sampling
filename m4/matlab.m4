
dnl  ====================================================
dnl             MATLAB library macros
dnl  ====================================================

dnl
dnl  CHECK_MATLAB - find Matlab libraries
dnl  ----------------------------------------------------
AC_DEFUN([CHECK_MATLAB],
[
AC_ARG_WITH(matlab,
	[  --with-matlab=DIR	the directory where Matlab is installed ],
	[if test -n "${withval}" -a ! "${withval}" = "yes"; then
		MATLAB_DIR=${withval};
	 else
		AC_PATH_PROG([MATLAB_PROG], [matlab])
		if test -n "${MATLAB_PROG}"; then
		   	MATLAB_DIR=$($MATLAB_PROG -e | grep '^MATLAB=' | sed s,MATLAB=,,);
		else
			AC_MSG_ERROR([You must specify the matlab directory with --with-matlab=PATH])
		fi;
	 fi],
	MATLAB_DIR=)

if test -n "${MATLAB_DIR}"
then
    AC_MSG_CHECKING(for Matlab software)

    case $build_os in
    *linux*)
	case $host_cpu in
	     ia64*|x86_64)
	     	MATLAB_FLAGS="-I\"${MATLAB_DIR}/extern/include\" -I\"${MATLAB_DIR}/simulink/include\" -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG";
       		MATLAB_LINK="-pthread -shared -Wl,--version-script,\"${MATLAB_DIR}/extern/lib/glnxa64/mexFunction.map\"";
		MATLAB_LINK_EXTRA="-Wl,--no-undefined ../src/.libs/libadaptive-sampling.a -static-libgcc";
       		MATLAB_LIB="-Wl,--rpath-link,\"${MATLAB_DIR}/extern/lib/glnxa64\",--rpath-link,\"${MATLAB_DIR}/bin/glnxa64\" -L\"${MATLAB_DIR}/bin/glnxa64\" -lmx -lmex -lmat -lm";
        	MEXEXT=mexa64;;
	     *)
	     	MATLAB_FLAGS="-I\"${MATLAB_DIR}/extern/include\" -I\"${MATLAB_DIR}/simulink/include\" -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG";
       		MATLAB_LINK="-pthread -shared -Wl,--version-script,\"${MATLAB_DIR}/extern/lib/glnx86/mexFunction.map\"";
		MATLAB_LINK_EXTRA="-Wl,--no-undefined ../src/.libs/libadaptive-sampling.a -static-libgcc";
       		MATLAB_LIB="-Wl,--rpath-link,\"${MATLAB_DIR}/extern/lib/glnx86\",--rpath-link,\"${MATLAB_DIR}/bin/glnx86\" -L\"${MATLAB_DIR}/bin/glnx86\" -lmx -lmex -lmat -lm";
        	MEXEXT=mexglx;;
	esac
	;;
    *darwin*)
	case $host_cpu in
	     ia64*|x86_64)
	     	MATLAB_FLAGS="-I\"${MATLAB_DIR}/extern/include\" -I\"${MATLAB_DIR}/simulink/include\" -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG";
       		MATLAB_LINK="-pthread -bundle -Wl,-flat_namespace -undefined suppress -Wl,-exported_symbols_list,\"${MATLAB_DIR}/extern/lib/maci64/mexFunction.map\"";
		MATLAB_LINK_EXTRA="../src/.libs/libadaptive-sampling.a";
       		MATLAB_LIB="-L\"${MATLAB_DIR}/bin/maci64\" -lmx -lmex -lmat -lm";
        	MEXEXT=mexmaci64;;
	     *powerpc*)
	     	MATLAB_FLAGS="-I\"${MATLAB_DIR}/extern/include\" -I\"${MATLAB_DIR}/simulink/include\" -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG";
       		MATLAB_LINK="-pthread -bundle -Wl,-flat_namespace -undefined suppress -Wl,-exported_symbols_list,\"${MATLAB_DIR}/extern/lib/mac/mexFunction.map\"";
		MATLAB_LINK_EXTRA="../src/.libs/libadaptive-sampling.a";
       		MATLAB_LIB="-L\"${MATLAB_DIR}/bin/mac\" -lmx -lmex -lmat -lm";
        	MEXEXT=mexmac;;
	     *)
	     	MATLAB_FLAGS="-I\"${MATLAB_DIR}/extern/include\" -I\"${MATLAB_DIR}/simulink/include\" -DMATLAB_MEX_FILE -fPIC -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG";
       		MATLAB_LINK="-pthread -bundle -Wl,-flat_namespace -undefined suppress -Wl,-exported_symbols_list,\"${MATLAB_DIR}/extern/lib/maci/mexFunction.map\"";
		MATLAB_LINK_EXTRA="../src/.libs/libadaptive-sampling.a";
       		MATLAB_LIB="-L\"${MATLAB_DIR}/bin/maci\" -lmx -lmex -lmat -lm";
        	MEXEXT=mexmaci;;
	esac
	;;
    *mingw*)
        MATLAB_FLAGS="-I\"${MATLAB_DIR}/extern/include\" -I\"${MATLAB_DIR}/simulink/include\" -fno-exceptions -DMATLAB_MEX_FILE -DNDEBUG";
        MATLAB_LINK="-shared -W1,--version-script,\"${MATLAB_DIR}/extern/lib/win32/mexFunction.def\"";
        MATLAB_LIB="-W1,--rpath-link,\"${MATLAB_DIR}/extern/lib/win32\",--rpath-link,\"${MATLAB_DIR}/bin/win32\" \"${MATLAB_DIR}/bin/win32/libmx.a\" \"${MATLAB_DIR}/bin/win32/libmex.a\" \"${MATLAB_DIR}/bin/win32/libmat.a\" -lm";
        MATLAB_LINK="-shared -L\"${MATLAB_DIR}/bin/win32\" -W1,--version-script,\"${MATLAB_DIR}/extern/lib/win32/mexFunction.def\"";
	MATLAB_LINK_EXTRA="-Wl,--no-undefined ../src/.libs/libadaptive-sampling.a -static-libgcc";
        MATLAB_LIB="-lmx -lmex -lmat -lm";
        MEXEXT=mexw32;
    esac
    AC_MSG_RESULT($MATLAB_LINK $MATLAB_LIB)
    AC_SUBST(MATLAB_DIR)
    AC_SUBST(MATLAB_LIB)
    AC_SUBST(MATLAB_LINK)
    AC_SUBST(MATLAB_LINK_EXTRA)
    AC_SUBST(MATLAB_FLAGS)
    AC_SUBST(MEXEXT)
    AC_SUBST(MEXVERSION)
    AM_CONDITIONAL(HAVE_MATLAB, test "xyes" = "xyes" )
    have_matlab=yes
else
    AM_CONDITIONAL(HAVE_MATLAB, test "xno" = "xyes" )
    have_matlab=no
fi
])
