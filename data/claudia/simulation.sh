#! /bin/sh

set -x

run_sampling() {
    adaptive-sampling -m -n 250 --save=simulation1-01.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-02.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-03.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-04.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-05.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-06.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-07.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-08.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-09.result simulation1.cfg
    adaptive-sampling -m -n 250 --save=simulation1-10.result simulation1.cfg

    adaptive-sampling -m -n 250 --save=simulation2-01.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-02.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-03.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-04.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-05.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-06.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-07.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-08.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-09.result simulation2.cfg
    adaptive-sampling -m -n 250 --save=simulation2-10.result simulation2.cfg

    adaptive-sampling -m -n 250 --save=simulation3-01.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-02.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-03.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-04.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-05.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-06.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-07.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-08.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-09.result simulation3.cfg
    adaptive-sampling -m -n 250 --save=simulation3-10.result simulation3.cfg

    adaptive-sampling -m -n 250 --save=simulation4-01.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-02.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-03.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-04.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-05.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-06.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-07.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-08.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-09.result simulation4.cfg
    adaptive-sampling -m -n 250 --save=simulation4-10.result simulation4.cfg
}

generate_plots() {
    adaptive-sampling -m --savefig=simulation1-prior.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-01.result --savefig=simulation1-01.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-02.result --savefig=simulation1-02.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-03.result --savefig=simulation1-03.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-04.result --savefig=simulation1-04.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-05.result --savefig=simulation1-05.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-06.result --savefig=simulation1-06.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-07.result --savefig=simulation1-07.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-08.result --savefig=simulation1-08.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-09.result --savefig=simulation1-09.pdf simulation1.cfg
    adaptive-sampling -m --load=simulation1-10.result --savefig=simulation1-10.pdf simulation1.cfg

    pdftk simulation1-prior.pdf simulation1-??.pdf cat output simulation1.pdf

    adaptive-sampling -m --savefig=simulation2-prior.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-01.result --savefig=simulation2-01.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-02.result --savefig=simulation2-02.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-03.result --savefig=simulation2-03.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-04.result --savefig=simulation2-04.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-05.result --savefig=simulation2-05.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-06.result --savefig=simulation2-06.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-07.result --savefig=simulation2-07.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-08.result --savefig=simulation2-08.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-09.result --savefig=simulation2-09.pdf simulation2.cfg
    adaptive-sampling -m --load=simulation2-10.result --savefig=simulation2-10.pdf simulation2.cfg

    pdftk simulation2-prior.pdf simulation2-??.pdf cat output simulation2.pdf

    adaptive-sampling -m --savefig=simulation3-prior.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-01.result --savefig=simulation3-01.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-02.result --savefig=simulation3-02.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-03.result --savefig=simulation3-03.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-04.result --savefig=simulation3-04.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-05.result --savefig=simulation3-05.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-06.result --savefig=simulation3-06.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-07.result --savefig=simulation3-07.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-08.result --savefig=simulation3-08.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-09.result --savefig=simulation3-09.pdf simulation3.cfg
    adaptive-sampling -m --load=simulation3-10.result --savefig=simulation3-10.pdf simulation3.cfg

    pdftk simulation3-prior.pdf simulation3-??.pdf cat output simulation3.pdf

    adaptive-sampling -m --savefig=simulation4-prior.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-01.result --savefig=simulation4-01.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-02.result --savefig=simulation4-02.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-03.result --savefig=simulation4-03.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-04.result --savefig=simulation4-04.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-05.result --savefig=simulation4-05.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-06.result --savefig=simulation4-06.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-07.result --savefig=simulation4-07.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-08.result --savefig=simulation4-08.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-09.result --savefig=simulation4-09.pdf simulation4.cfg
    adaptive-sampling -m --load=simulation4-10.result --savefig=simulation4-10.pdf simulation4.cfg

    pdftk simulation4-prior.pdf simulation4-??.pdf cat output simulation4.pdf
}

case "$1" in
    sample)
        run_sampling
        ;;
    plot)
        generate_plots
        ;;
esac
