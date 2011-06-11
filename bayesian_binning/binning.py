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
    import bayesian_binning.visualization as vis

import bayesian_binning.config     as config
import bayesian_binning.interface  as interface
import bayesian_binning.statistics as statistics

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
    print "   -b                             - compute break probabilities"
    print "   -d                             - compute differential gain"
    print "       --effective-counts         - compute effective counts"
    print "   -m  --marginal                 - compute full marginal distribution"
    print "   -r  --marginal-range=(FROM,TO) - limit range for the marginal distribution"
    print "   -s  --marginal-step=STEP       - step size for the marginal distribution"
    print "       --epsilon=EPSILON          - epsilon for entropy estimations"
    print "   -k  --moments=N                - compute the first N>=3 moments"
    print "       --which=EVENT              - for which event to compute the binning"
    print
    print "       --threads=THREADS          - number of threads [default: 1]"
    print
    print "       --load=FILE                - load result from file"
    print "       --save=FILE                - save result to file"
    print "       --savefig=FILE             - save figure to file"
    print
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print "   -t, --prombsTest               - test prombs algorithm"
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
        'mpost'     : mpost,
        'differential_gain' : [] }
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
        options['visualization'] = config.readVisualization(config_parser, 'Counts', os.path.dirname(config_file))
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
        options['visualization'] = config.readVisualization(config_parser, 'Trials', os.path.dirname(config_file))
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
    'epsilon'           : 0.00001,
    'marginal'          : 0,
    'marginal_step'     : 0.01,
    'marginal_range'    : (0.0,1.0),
    'n_moments'         : 3,
    'which'             : 0,
    'threads'           : 1,
    'visualization'     : None,
    'load'              : None,
    'save'              : None,
    'savefig'           : None,
    'verbose'           : False,
    'prombsTest'        : False,
    'compare'           : False,
    'bprob'             : False,
    'differential_gain' : False,
    'effective_counts'  : False,
    'model_posterior'   : True,
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=", "marginal", "marginal-range:"
                      "marginal-step=", "which=", "epsilon=", "moments=", "prombsTest",
                      "effective-counts", "savefig=", "threads="]
        opts, tail = getopt.getopt(sys.argv[1:], "dmr:s:k:bhvt", longopts)
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
            if int(a) >= 3:
                options["n_moments"] = int(a)
            else:
                usage()
                return 0
        if o in ("-h", "--help"):
            usage()
            return 0
        if o == "-b":
            options["bprob"] = True
        if o == "-d":
            options["differential_gain"] = True
        if o == "--effective-counts":
            options["effective_counts"] = True
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
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
