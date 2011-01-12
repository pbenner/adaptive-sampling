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
    print "       --epsilon=EPSILON       - epsilon for entropy estimations"
    print "   -m  --moments=N             - compute the first N>=3 moments"
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
    config = ConfigParser.RawConfigParser()
    config.read(options['load'])
    if not config.has_section('Result'):
        raise IOError("Invalid configuration file.")

    moments_str = config.get   ('Result', 'moments')
    moments     = []
    for line in moments_str.split('\n'):
        if line != '':
            moments.append([float(a) for a in line.split(' ')])
    mpost_str   = config.get('Result', 'mpost')
    mpost       = map(float, mpost_str.split(' '))
    if config.has_option('Result', 'bprob'):
        bprob_str = config.get('Result', 'bprob')
        bprob     = map(float, bprob_str.split(' '))
    else:
        bprob     = []
    if config.has_option('Result', 'multibin_entropy'):
        multibin_entropy_str = config.get('Result', 'multibin_entropy')
        multibin_entropy     = map(float, multibin_entropy_str.split(' '))
    else:
        multibin_entropy     = []
    result = {
        'moments' : moments,
        'bprob'   : bprob,
        'mpost'   : mpost,
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
    config.set('Result', 'bprob',   " ".join(map(str, result['bprob'])))
    config.set('Result', 'mpost',   " ".join(map(str, result['mpost'])))
    config.set('Result', 'multibin_entropy', " ".join(map(str, result['multibin_entropy'])))
    configfile = open(options['save'], 'wb')
    config.write(configfile)

# parse config file
# ------------------------------------------------------------------------------

def timingsToCounts(timings, binsize):
    MIN    = min(map(min, timings))
    MAX    = max(map(max, timings))
    N      = int(math.ceil(float(MAX-MIN)/binsize))
    counts = list(np.repeat(0, N+1))
    x      = range(MIN, MAX+binsize, binsize)
    for trial in timings:
        for t in trial:
            n = int(math.ceil(float(t-MIN)/binsize))
            counts[n] += 1
    return x, counts

def computeFailures(successes, trials):
    N        = len(successes)
    failures = np.repeat(trials, N) - successes
    if any([ a<0 for a in failures]):
        raise ValueError("Number of trials is smaller than some counts.")
    return failures

def readModelPrior(models_str, N):
    models = []
    mprior = list(np.repeat(0, N))
    for str in models_str.split(' '):
        models.append(int(str))
    num_models = len(models)
    for model in models:
        mprior[model-1] = 1.0/num_models
    return mprior

def readAlpha(config, section, n):
    alpha = []
    if config.has_option(section, 'alpha'):
        alpha_str  = config.get   (section, 'alpha')
        for str in alpha_str.split(' '):
            alpha.append(int(str))
        if len(alpha) < n:
            raise ValueError("Not enough alpha parameters.")
        if len(alpha) > n:
            raise ValueError("Too many alpha parameters.")
    else:
        alpha = list(np.repeat(1, n))
    return alpha

def parseConfig(file):
    config = ConfigParser.RawConfigParser()
    config.read(file)

    if config.sections() == []:
        raise IOError("Invalid configuration file.")
    if config.has_section('Counts'):
        counts_str  = config.get   ('Counts', 'counts')
        counts      = []
        for line in counts_str.split('\n'):
            if line != '':
                counts.append([int(a) for a in line.split(' ')])
        N           = len(counts[0])
        prior       = list(np.repeat(1, N))
        alpha       = readAlpha(config, 'Counts', len(counts))
        if config.has_option('Counts', 'mprior'):
            prior   = readModelPrior(config.get('Counts', 'mprior'), N)
        result      = bin(counts, alpha, prior)
        if options['save']:
            saveResult(result)
        else:
            x = np.arange(0, N, 1)
            vis.plotBinning(x, result, options['bprob'], options['multibin_entropy'])
    if config.has_section('Trials'):
        binsize     = config.getint('Trials', 'binsize')
        timings_str = config.get   ('Trials', 'timings')
        timings     = []
        for line in timings_str.split('\n'):
            if line != '':
                timings.append([int(a) for a in line.split(' ')])

        x, successes = timingsToCounts(timings, binsize)
        trials       = len(timings)
        failures     = computeFailures(successes, trials)
        counts       = [successes, failures]
        N            = len(successes)
        prior        = list(np.repeat(1, N))
        alpha        = readAlpha(config, 'Trials', 2)
        if config.has_option('Trials', 'mprior'):
            prior    = readModelPrior(config.get('Trials', 'mprior'), N)
        result       = bin(counts, alpha, prior)
        if options['save']:
            saveResult(result)
        else:
            vis.plotBinningSpikes(x, timings, result, options['bprob'], options['multibin_entropy'])

# main
# ------------------------------------------------------------------------------

options = {
    'epsilon'    : 0.00001,
    'n_moments'  : 3,
    'which'      : 0,
    'load'       : None,
    'save'       : None,
    'verbose'    : False,
    'prombsTest' : False,
    'compare'    : False,
    'bprob'      : False,
    'differential_gain' : False,
    'multibin_entropy'  : False,
    'model_posterior'   : True,
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=",
                      "which=", "epsilon=", "moments=", "prombsTest"]
        opts, tail = getopt.getopt(sys.argv[1:], "dem:bhvt", longopts)
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
        if o in ("-m", "--moments"):
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
