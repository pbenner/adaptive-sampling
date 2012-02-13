#! /bin/sh

if [ ! $# -eq 1 ]; then
    echo "Usage: matlab.sh MATLAB_DIR"
    exit 1
fi

MATLAB=$1
MAPFILE=mexFunction.map
ARCH_LIST='glnx86 glnxa64 mac maci maci64 sol2 sol64'

# check_archlist 
# ------------------------------------------------------------------------------
check_archlist () {
    if [ $# -gt 0 ]; then
        arch_in=`expr "$1" : '.*=\(.*\)'`
        if [ "$arch_in" != "" ]; then
            ARCH=`echo "$ARCH_LIST EOF $arch_in" | awk '
# ------------------------------------------------------------------------------
        { for (i = 1; i <= NF; i = i + 1)
              if ($i == "EOF")
                  narch = i - 1
          for (i = 1; i <= narch; i = i + 1)
                if ($i == $NF || "-" $i == $NF) {
                    print $i
                    exit
                }
        }'`
# ------------------------------------------------------------------------------
            if [ "$ARCH" = "" -a $# -eq 1 ]; then
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
echo ' '
echo "^G    Warning: $1 does not specify a valid architecture - ignored . . ."
echo ' '
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            fi
        else
            ARCH=""
        fi
    else
        ARCH=""
    fi
#
    return 0
}
# ------------------------------------------------------------------------------

. $MATLAB/bin/util/arch.sh
. $MATLAB/bin/mexopts.sh

# remove -pthread since clang doesn't support it
CFLAGS=$(echo $CFLAGS | sed -e s/-pthread//)
LDFLAGS=$(echo $LDFLAGS | sed -e s/-pthread//)

echo MATLAB_RPATH=\'$RPATH\'\;
echo MATLAB_CFLAGS=\'-I\"$MATLAB/extern/include\" $CFLAGS\'\;
echo MATLAB_LIBS=\'$CLIBS\'\;
echo MATLAB_LDFLAGS=\'$LDFLAGS\'\;
echo MATLAB_LDEXTENSION=\'$LDEXTENSION\'\;
