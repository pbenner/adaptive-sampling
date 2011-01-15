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
from matplotlib import *
from matplotlib.pyplot import *
import matplotlib.patches as patches
import matplotlib.path as path

import bayesian_binning.config as config
import bayesian_binning.interface as interface
import bayesian_binning.visualization as vis
import bayesian_binning.statistics as statistics

# global options
# ------------------------------------------------------------------------------

verbose = False

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "bayesian-binning.py [option]... FILE "
    print
    print "Options:"
    print "   -b                          - compute break probabilities"
    print "   -d                          - compute differential gain"
    print "   -e                          - compute multibin entropies"
    print "   -m  --marginal              - compute full marginal distribution"
    print "   -s  --marginal-step=STEP    - step size for the marginal distribution"
    print "       --epsilon=EPSILON       - epsilon for entropy estimations"
    print "   -k  --moments=N             - compute the first N>=3 moments"
    print "       --which=EVENT           - for which event to compute the binning"
    print
    print "       --load                  - load result from file"
    print "       --save                  - save result to file"
    print
    print "   -h, --help                  - print help"
    print "   -v, --verbose               - be verbose"
    print "   -t, --prombsTest            - test prombs algorithm"
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
    if config_parser.has_option('Result', 'multibin_entropy'):
        multibin_entropy = config.readVector(config_parser, 'Result', 'multibin_entropy', float)
    else:
        multibin_entropy = []
    result = {
        'moments'   : moments,
        'marginals' : marginals,
        'bprob'     : bprob,
        'mpost'     : mpost,
        'multibin_entropy'  : multibin_entropy,
        'differential_gain' : [] }
    return result

# binning
# ------------------------------------------------------------------------------

def bin(counts, alpha, mprior):
    """Call the binning library."""
    if options['load']:
        return load_config()
    else:
        counts_i = [ map(int, row) for row in counts ]
        mprior_i =   map(float, mprior)
        alpha_i  =   map(int, alpha)
        return interface.binning(counts_i, alpha_i, mprior_i, options)

# save result
# ------------------------------------------------------------------------------

def saveResult(result):
    config = ConfigParser.ConfigParser()
    config.add_section('Result')
    config.set('Result', 'moments', "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['moments'])))
    config.set('Result', 'marginals', "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['marginals'])))
    config.set('Result', 'bprob',   " ".join(map(str, result['bprob'])))
    config.set('Result', 'mpost',   " ".join(map(str, result['mpost'])))
    config.set('Result', 'multibin_entropy', " ".join(map(str, result['multibin_entropy'])))
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

def timingsToCounts(timings, binsize):
    MIN       = min(map(min, timings))
    MAX       = max(map(max, timings))
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
    return x, counts

def parseConfig(config_file):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)

    if config_parser.sections() == []:
        raise IOError("Invalid configuration file.")
    if config_parser.has_section('Counts'):
        counts    = config.readMatrix(config_parser, 'Counts', 'counts', int)
        N         = len(counts[0])
        prior     = list(np.repeat(1, N))
        alpha     = config.readAlpha(config_parser, len(counts), 'Counts', int)
        prior     = config.readModelPrior(config_parser, N, 'Counts', int)
        result    = bin(counts, alpha, prior)
        if options['save']:
            saveResult(result)
        else:
            x = np.arange(0, N, 1)
            vis.plotBinning(x, result, options)
    if config_parser.has_section('Trials'):
        binsize   = config_parser.getint('Trials', 'binsize')
        timings   = config.readMatrix(config_parser, 'Trials', 'timings', int)
        x, counts = timingsToCounts(timings, binsize)
        N         = len(counts[0])
        prior     = list(np.repeat(1, N))
        alpha     = config.readAlpha(config_parser, len(counts), 'Trials', int)
        prior     = config.readModelPrior(config_parser, N, 'Trials', int)
        result    = bin(counts, alpha, prior)
        if options['save']:
            saveResult(result)
        else:
            vis.plotBinningSpikes(x, timings, result, options)

# main
# ------------------------------------------------------------------------------

options = {
    'epsilon'       : 0.00001,
    'marginal'      : 0,
    'marginal_step' : 0.01,
    'n_moments'     : 3,
    'which'         : 0,
    'load'          : None,
    'save'          : None,
    'verbose'       : False,
    'prombsTest'    : False,
    'compare'       : False,
    'bprob'         : False,
    'differential_gain' : False,
    'multibin_entropy'  : False,
    'model_posterior'   : True,
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=", "marginal",
                      "marginal-step=", "which=", "epsilon=", "moments=", "prombsTest"]
        opts, tail = getopt.getopt(sys.argv[1:], "dems:k:bhvt", longopts)
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
        if o == "-e":
            options["multibin_entropy"] = True
        if o == "--load":
            options["load"] = a
        if o == "--save":
            options["save"] = a
        if o == "--which":
            options["which"] = int(a)
        if o == "--epsilon":
            options["epsilon"] = float(a)
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
