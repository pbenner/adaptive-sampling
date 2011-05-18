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

    adaptive-sampling -m -n 250 --save=simulation5-01.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-02.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-03.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-04.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-05.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-06.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-07.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-08.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-09.result simulation5.cfg
    adaptive-sampling -m -n 250 --save=simulation5-10.result simulation5.cfg

    adaptive-sampling -m -n 250 --save=simulation6-01.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-02.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-03.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-04.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-05.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-06.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-07.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-08.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-09.result simulation6.cfg
    adaptive-sampling -m -n 250 --save=simulation6-10.result simulation6.cfg

    adaptive-sampling -m -n 60  --strategy=uniform --save=simulation7-00.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-01.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-02.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-03.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-04.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-05.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-06.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-07.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-08.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-09.result simulation7.cfg
    adaptive-sampling -m -n 190 --strategy=differential-gain --load=simulation7-00.result --save=simulation7-10.result simulation7.cfg

    adaptive-sampling -m -n 250 --save=simulation8-01.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-02.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-03.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-04.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-05.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-06.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-07.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-08.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-09.result simulation8.cfg
    adaptive-sampling -m -n 250 --save=simulation8-10.result simulation8.cfg

    adaptive-sampling -m -n 250 --save=simulation9-01.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-02.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-03.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-04.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-05.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-06.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-07.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-08.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-09.result simulation9.cfg
    adaptive-sampling -m -n 250 --save=simulation9-10.result simulation9.cfg

    adaptive-sampling -m -n 250 --save=simulation10-01.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-02.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-03.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-04.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-05.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-06.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-07.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-08.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-09.result simulation10.cfg
    adaptive-sampling -m -n 250 --save=simulation10-10.result simulation10.cfg
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

    adaptive-sampling -m --savefig=simulation5-prior.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-01.result --savefig=simulation5-01.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-02.result --savefig=simulation5-02.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-03.result --savefig=simulation5-03.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-04.result --savefig=simulation5-04.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-05.result --savefig=simulation5-05.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-06.result --savefig=simulation5-06.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-07.result --savefig=simulation5-07.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-08.result --savefig=simulation5-08.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-09.result --savefig=simulation5-09.pdf simulation5.cfg
    adaptive-sampling -m --load=simulation5-10.result --savefig=simulation5-10.pdf simulation5.cfg

    pdftk simulation5-prior.pdf simulation5-??.pdf cat output simulation5.pdf

    adaptive-sampling -m --savefig=simulation6-prior.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-01.result --savefig=simulation6-01.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-02.result --savefig=simulation6-02.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-03.result --savefig=simulation6-03.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-04.result --savefig=simulation6-04.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-05.result --savefig=simulation6-05.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-06.result --savefig=simulation6-06.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-07.result --savefig=simulation6-07.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-08.result --savefig=simulation6-08.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-09.result --savefig=simulation6-09.pdf simulation6.cfg
    adaptive-sampling -m --load=simulation6-10.result --savefig=simulation6-10.pdf simulation6.cfg

    pdftk simulation6-prior.pdf simulation6-??.pdf cat output simulation6.pdf

    adaptive-sampling -m --savefig=simulation7-prior.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-01.result --savefig=simulation7-01.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-02.result --savefig=simulation7-02.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-03.result --savefig=simulation7-03.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-04.result --savefig=simulation7-04.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-05.result --savefig=simulation7-05.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-06.result --savefig=simulation7-06.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-07.result --savefig=simulation7-07.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-08.result --savefig=simulation7-08.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-09.result --savefig=simulation7-09.pdf simulation7.cfg
    adaptive-sampling -m --load=simulation7-10.result --savefig=simulation7-10.pdf simulation7.cfg

    pdftk simulation7-prior.pdf simulation7-??.pdf cat output simulation7.pdf

    adaptive-sampling -m --savefig=simulation8-prior.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-01.result --savefig=simulation8-01.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-02.result --savefig=simulation8-02.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-03.result --savefig=simulation8-03.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-04.result --savefig=simulation8-04.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-05.result --savefig=simulation8-05.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-06.result --savefig=simulation8-06.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-07.result --savefig=simulation8-07.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-08.result --savefig=simulation8-08.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-09.result --savefig=simulation8-09.pdf simulation8.cfg
    adaptive-sampling -m --load=simulation8-10.result --savefig=simulation8-10.pdf simulation8.cfg

    pdftk simulation8-prior.pdf simulation8-??.pdf cat output simulation8.pdf

    adaptive-sampling -m --savefig=simulation9-prior.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-01.result --savefig=simulation9-01.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-02.result --savefig=simulation9-02.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-03.result --savefig=simulation9-03.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-04.result --savefig=simulation9-04.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-05.result --savefig=simulation9-05.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-06.result --savefig=simulation9-06.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-07.result --savefig=simulation9-07.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-08.result --savefig=simulation9-08.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-09.result --savefig=simulation9-09.pdf simulation9.cfg
    adaptive-sampling -m --load=simulation9-10.result --savefig=simulation9-10.pdf simulation9.cfg

    pdftk simulation9-prior.pdf simulation9-??.pdf cat output simulation9.pdf

    adaptive-sampling -m --savefig=simulation10-prior.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-01.result --savefig=simulation10-01.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-02.result --savefig=simulation10-02.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-03.result --savefig=simulation10-03.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-04.result --savefig=simulation10-04.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-05.result --savefig=simulation10-05.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-06.result --savefig=simulation10-06.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-07.result --savefig=simulation10-07.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-08.result --savefig=simulation10-08.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-09.result --savefig=simulation10-09.pdf simulation10.cfg
    adaptive-sampling -m --load=simulation10-10.result --savefig=simulation10-10.pdf simulation10.cfg

    pdftk simulation10-prior.pdf simulation10-??.pdf cat output simulation10.pdf
}

case "$1" in
    sample)
        run_sampling
        ;;
    plot)
        generate_plots
        ;;
esac
