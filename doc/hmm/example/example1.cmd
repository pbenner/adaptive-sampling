#! /bin/sh

set -x

run_sampling() {
    adaptive-sampling -n   1 -m --hmm --rho 0.8 --save=example1.result1 example1.cfg
    adaptive-sampling -n   1 -m --hmm --rho 0.8 --load=example1.result1 --save=example1.result2 example1.cfg
    adaptive-sampling -n   8 -m --hmm --rho 0.8 --load=example1.result2 --save=example1.result3 example1.cfg
    adaptive-sampling -n  40 -m --hmm --rho 0.8 --load=example1.result3 --save=example1.result4 example1.cfg
    adaptive-sampling -n  50 -m --hmm --rho 0.8 --load=example1.result4 --save=example1.result5 example1.cfg
    adaptive-sampling -n 100 -m --hmm --rho 0.8 --load=example1.result5 --save=example1.result6 example1.cfg
}

generate_plots() {
    adaptive-sampling -m --hmm --rho 0.8 --no-model-posterior --savefig=example1-0.pdf example1.cfg
    adaptive-sampling -m --hmm --rho 0.8 --no-model-posterior --load=example1.result1 --savefig=example1-1.pdf example1.cfg
    adaptive-sampling -m --hmm --rho 0.8 --no-model-posterior --load=example1.result2 --savefig=example1-2.pdf example1.cfg
    adaptive-sampling -m --hmm --rho 0.8 --no-model-posterior --load=example1.result3 --savefig=example1-3.pdf example1.cfg
    adaptive-sampling -m --hmm --rho 0.8 --no-model-posterior --load=example1.result4 --savefig=example1-4.pdf example1.cfg
    adaptive-sampling -m --hmm --rho 0.8 --no-model-posterior --load=example1.result5 --savefig=example1-5.pdf example1.cfg
    adaptive-sampling -m --hmm --rho 0.8 --no-model-posterior --load=example1.result6 --savefig=example1-6.pdf example1.cfg

    pdftk example1-?.pdf cat output example1.pdf
    rm -f example1-?.pdf
}

case "$1" in
    sample)
	run_sampling
	;;
    plot)
	generate_plots
	;;
esac
