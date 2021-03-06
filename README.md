## Installation

For local installations, say $HOME/.usr, remember to set the following
environment variables:

	export PATH="$HOME/.usr/bin:$PATH"
	export CPATH=$HOME/.usr/include
	export LD_LIBRARY_PATH=$HOME/.usr/lib:$LD_LIBRARY_PATH
	export LIBRARY_PATH=$LD_LIBRARY_PATH
	export MANPATH=$HOME/.usr/share/man:$MANPATH

and then compile with

	./configure --prefix=$HOME/.usr
	make
	make install


## Binning examples

	bayesian-binning -v -b -m data/binning/data1.cfg
	bayesian-binning -m --which=1 data/binning/data1.cfg

	bayesian-binning -v -b -m data/binning-spikes/data1.cfg

	bayesian-binning -v -b -m --save=data/binning-spikes/data1.result data/binning-spikes/data1.cfg
	bayesian-binning -v -b -m --load=data/binning-spikes/data1.result data/binning-spikes/data1.cfg

	bayesian-binning -v -b -m -r '(0,0.1)' -s 0.001 --save=data/binning-spikes/data1.result data/binning-spikes/data1.cfg

	bayesian-binning -v -b --threads=42 G833C11R4_2tone-PSTH_Exci.cfg


## Adaptive sampling examples

	adaptive-sampling -n  10 -m data/sampling/data1.cfg

	adaptive-sampling -n  10 -m --save=data/sampling/data1.result1 data/sampling/data1.cfg
	adaptive-sampling -n  90 -m --load=data/sampling/data1.result1 --save=data/sampling/data1.result2 data/sampling/data1.cfg
	adaptive-sampling -n 100 -m --load=data/sampling/data1.result2 --save=data/sampling/data1.result3 data/sampling/data1.cfg

	adaptive-sampling -m --load=data/sampling/data1.result3 data/sampling/data1.cfg

	adaptive-sampling -n 100 -m --plot-utility data/sampling/data1.cfg

	adaptive-sampling -m -n 400 --strategy=uniform --save=data/sampling/data1.result4 data/sampling/data1.cfg

	adaptive-sampling -n 5 -m --strategy=effective-counts data/sampling/data4.cfg 2>/dev/null
	adaptive-sampling -n 5 -m --strategy=kl-divergence data/sampling/data4.cfg 2>/dev/null
	adaptive-sampling -n 5 -m --strategy=kl-divergence --kl-multibin data/sampling/data4.cfg 2>/dev/null
