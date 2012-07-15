#! /usr/bin/env python

# Copyright (C) 2010 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys
import getopt
import os
import ConfigParser
import numpy as np
import math

from itertools import izip

def importMatplotlib(backend=None):
    global vis
    from matplotlib import use
    if backend:
        use(backend)
    import visualization as vis

import config
import interface
import statistics

# global options
# ------------------------------------------------------------------------------

verbose = False

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "bayesian-binning [option]... FILE "
    print
    print "Options:"
    print "   -b                                - compute break probabilities"
    print "       --hmm                         - use hidden Markov model"
    print "   -m  --marginal                    - compute full marginal distribution"
    print "   -r  --marginal-range=(FROM,TO)    - limit range for the marginal distribution"
    print "   -s  --marginal-step=STEP          - step size for the marginal distribution"
    print "       --no-model-posterior          - do not compute the model posterior"
    print "       --epsilon=EPSILON             - epsilon for the extended prombs"
    print "   -k  --moments=N                   - compute the first N>=2 moments"
    print "       --which=EVENT                 - for which event to compute the binning"
    print "       --algorithm=NAME              - select an algorithm [mgs, prombstree, default: prombs]"
    print "       --mgs-samples=BURN_IN:SAMPLES - number of samples [default: 100:2000]"
    print
    print "       --threads=THREADS             - number of threads [default: 1]"
    print "       --stacksize=BYTES             - thread stack size [default: 256*1024]"
    print
    print "       --load=FILE                   - load result from file"
    print "       --save=FILE                   - save result to file"
    print "       --savefig=FILE                - save figure to file"
    print
    print "   -h, --help                        - print help"
    print "   -v, --verbose                     - be verbose"
    print "   -t, --prombsTest                  - test prombs algorithm"
    print

# load results from file
# ------------------------------------------------------------------------------

def load_config():
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(options['load'])
    if not config_parser.has_section('Result'):
        raise IOError("Invalid configuration file.")

    moments   = config.readMatrix(config_parser, 'Result', 'moments',   float)
    marginals = config.readMatrix(config_parser, 'Result', 'marginals', float)
    mpost     = config.readVector(config_parser, 'Result', 'mpost',     float)
    if config_parser.has_option('Result', 'bprob'):
        bprob = config.readVector(config_parser, 'Result', 'bprob',     float)
    else:
        bprob     = []
    result = {
        'moments'   : moments,
        'marginals' : marginals,
        'bprob'     : bprob,
        'mpost'     : mpost }
    return result

# binning
# ------------------------------------------------------------------------------

def bin(counts, alpha, beta, gamma):
    """Call the binning library."""
    if options['load']:
        return load_config()
    else:
        events = len(counts)
        return interface.binning(events, counts, alpha, beta, gamma, options)

# save result
# ------------------------------------------------------------------------------

def saveResult(result):
    config = ConfigParser.ConfigParser()
    config.add_section('Result')
    config.set('Result', 'moments', "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['moments'])))
    config.set('Result', 'marginals', "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['marginals'])))
    config.set('Result', 'bprob',   " ".join(map(str, result['bprob'])))
    config.set('Result', 'mpost',   " ".join(map(str, result['mpost'])))
    configfile = open(options['save'], 'wb')
    config.write(configfile)

# parse config file
# ------------------------------------------------------------------------------

def computeFailures(successes, trials):
    N        = len(successes)
    failures = np.repeat(trials, N) - successes
    if any([ a<0 for a in failures]):
        raise ValueError("Number of trials is smaller than some counts.")
    return failures

def timingsToCounts(timings, binsize, srange):
    if srange == None:
        MIN   = min(map(min, timings))
        MAX   = max(map(max, timings))
    else:
        MIN   = srange[0]
        MAX   = srange[1]
    N         = int(math.ceil(float(MAX-MIN)/binsize))
    successes = list(np.repeat(0, N+1))
    x         = range(MIN, MAX+binsize, binsize)
    for trial in timings:
        for t in trial:
            n = int(math.ceil(float(t-MIN)/binsize))
            successes[n] += 1
    trials    = len(timings)
    failures  = computeFailures(successes, trials)
    counts    = [successes, failures]
    return x, statistics.countStatistic(counts)

def parseConfig(config_file):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)

    if config_parser.sections() == []:
        raise IOError("Invalid configuration file.")
    if config_parser.has_section('Counts'):
        config.readVisualization(config_parser, 'Counts', os.path.dirname(config_file), options)
        config.readAlgorithm(config_parser, 'Counts', os.path.dirname(config_file), options)
        config.readMgsSamples(config_parser, 'Counts', os.path.dirname(config_file), options)
        counts = config.readCounts(config_parser, 'Counts')
        K, L   = len(counts), len(counts[0])
        alpha, beta, gamma = config.getParameters(config_parser, 'Counts', os.path.dirname(config_file), K, L)
        result = bin(counts, alpha, beta, gamma)
        if options['save']:
            saveResult(result)
        else:
            if options['savefig']:
                importMatplotlib('Agg')
                from matplotlib.pyplot import savefig
                vis.plotBinning(result, options)
                savefig(options['savefig'], bbox_inches='tight', pad_inches=0)
            else:
                importMatplotlib()
                from matplotlib.pyplot import show
                vis.plotBinning(result, options)
                show()
    if config_parser.has_section('Trials'):
        config.readVisualization(config_parser, 'Trials', os.path.dirname(config_file), options)
        config.readAlgorithm(config_parser, 'Trials', os.path.dirname(config_file), options)
        config.readMgsSamples(config_parser, 'Trials', os.path.dirname(config_file), options)
        binsize   = config_parser.getint('Trials', 'binsize')
        timings   = config.readMatrix(config_parser, 'Trials', 'timings', int)
        srange    = None
        if config_parser.has_option('Trials', 'range'):
            srange = config.readVector(config_parser, 'Trials', 'range', int)
        x, counts = timingsToCounts(timings, binsize, srange)
        K, L   = len(counts), len(counts[0])
        alpha, beta, gamma = config.getParameters(config_parser, 'Trials', os.path.dirname(config_file), K, L)
        result    = bin(counts, alpha, beta, gamma)
        if options['save']:
            saveResult(result)
        else:
            if options['savefig']:
                importMatplotlib('Agg')
                from matplotlib.pyplot import savefig
                vis.plotBinningSpikes(x, timings, result, options)
                savefig(options['savefig'], bbox_inches='tight', pad_inches=0)
            else:
                importMatplotlib()
                from matplotlib.pyplot import show
                vis.plotBinningSpikes(x, timings, result, options)
                show()

# main
# ------------------------------------------------------------------------------

options = {
    'epsilon'              : 0.00001,
    'mgs_samples'          : (100,2000),
    'marginal'             : 0,
    'marginal_step'        : 0.01,
    'marginal_range'       : (0.0,1.0),
    'n_moments'            : 2,
    'which'                : 0,
    'threads'              : 1,
    'stacksize'            : 256*1024,
    'algorithm'            : 'prombs',
    'visualization'        : None,
    'load'                 : None,
    'save'                 : None,
    'savefig'              : None,
    'verbose'              : False,
    'prombsTest'           : False,
    'compare'              : False,
    'bprob'                : False,
    'utility'              : False,
    'kl_component'         : False,
    'kl_multibin'          : False,
    'effective_counts'     : False,
    'effective_posterior_counts' : False,
    'model_posterior'      : True,
    'hmm'                  : False,
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=", "marginal", "marginal-range:"
                      "marginal-step=", "which=", "epsilon=", "moments=", "prombsTest",
                      "savefig=", "threads=", "stacksize=", "algorithm=",
                      "mgs-samples=", "no-model-posterior", "hmm"]
        opts, tail = getopt.getopt(sys.argv[1:], "mr:s:k:bhvt", longopts)
    except getopt.GetoptError:
        usage()
        return 2
    output = None
    for o, a in opts:
        if o in ("-v", "--verbose"):
            sys.stderr.write("Verbose mode turned on.\n")
            options["verbose"] = True
        if o in ("-t", "--prombsTest"):
            sys.stderr.write("Testing prombs.\n")
            options["prombsTest"] = True
        if o in ("-m", "--marginal"):
            options["marginal"] = True
        if o in ("-r", "--marginal-range"):
            options['marginal_range'] = tuple(map(float, a[1:-1].split(',')))
        if o in ("-s", "--marginal-step"):
            options["marginal_step"] = float(a)
        if o in ("-k", "--moments"):
            if int(a) >= 2:
                options["n_moments"] = int(a)
            else:
                usage()
                return 0
        if o in ("-h", "--help"):
            usage()
            return 0
        if o == "-b":
            options["bprob"] = True
        if o == "--load":
            options["load"] = a
        if o == "--save":
            options["save"] = a
        if o == "--savefig":
            options["savefig"] = a
        if o == "--which":
            options["which"] = int(a)
        if o == "--epsilon":
            options["epsilon"] = float(a)
        if o == "--threads":
            if int(a) >= 1:
                options["threads"] = int(a)
            else:
                usage()
                return 0
        if o == "--stacksize":
            if int(a) >= 1:
                options["stacksize"] = int(a)
            else:
                usage()
                return 0
        if o == "--algorithm":
            options["algorithm"] = a
        if o == "--mgs-samples":
            options["mgs_samples"] = tuple(map(int, a.split(":")))
        if o == "--no-model-posterior":
            options["model_posterior"] = False
        if o == "--hmm":
            options["hmm"] = True
    if len(tail) != 1:
        usage()
        return 1
    interface.init(options['epsilon'])
    parseConfig(tail[0])
    interface.free()
    return 0

if __name__ == "__main__":
    sys.exit(main())
