#! /bin/sh

set -x

run_sampling() {
    echo
}

generate_plots() {
    ./example2.py

    pdftk example2-?.pdf cat output example2.pdf
}


case "$1" in
    sample)
	run_sampling
	;;
    plot)
	generate_plots
	;;
esac
