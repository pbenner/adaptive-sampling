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
import random
from itertools import izip
from matplotlib import *
from matplotlib.pyplot import *
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
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
    print "adaptive-sampling.py [option]... FILE "
    print
    print "Options:"
    print "   -b                             - compute break probabilities"
    print "       --blocks=N                 - sample in blocks of N measurements"
    print "   -d                             - compute differential gain"
    print "   -m  --marginal                 - compute full marginal distribution"
    print "   -r  --marginal-range=(FROM,TO) - limit range for the marginal distribution"
    print "   -s  --marginal-step=STEP       - step size for the marginal distribution"
    print "       --epsilon=EPSILON          - epsilon for entropy estimations"
    print "   -n  --samples=N                - number of samples"
    print "   -k  --moments=N                - compute the first N>=3 moments"
    print "       --strategy=STRATEGY        - uniform, uniform-random, differential-gain (default),"
    print "                                    effective-counts, or variance"
    print "       --which=EVENT              - for which event to compute the binning"
    print
    print "       --plot-utility             - plot utility as a function of sample steps"
    print
    print "       --load=FILE                - load result from file"
    print "       --save=FILE                - save result to file"
    print
    print "   -h, --help                     - print help"
    print "   -v, --verbose                  - be verbose"
    print "   -t, --prombsTest               - test prombs algorithm"
    print

# tools
# ------------------------------------------------------------------------------

def argmin(array):
    result = []
    for i in izip(array, xrange(len(array))):
        if i[0] == min(array):
            result.append(i[1])
    return result

def argmax(array):
    result = []
    for i in izip(array, xrange(len(array))):
        if i[0] == max(array):
            result.append(i[1])
    return result

def selectRandom(array):
    length = len(array)-1
    index  = random.randint(0, length)
    return array[index]

# save result
# ------------------------------------------------------------------------------

def saveResult(result):
    config = ConfigParser.ConfigParser()
    config.add_section('Sampling Result')
    config.set('Sampling Result', 'counts',    "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['counts'])))
    config.set('Sampling Result', 'moments',   "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['moments'])))
    config.set('Sampling Result', 'marginals', "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['marginals'])))
    config.set('Sampling Result', 'samples',   " ".join(map(str, result['samples'])))
    config.set('Sampling Result', 'bprob',     " ".join(map(str, result['bprob'])))
    config.set('Sampling Result', 'mpost',     " ".join(map(str, result['mpost'])))
    configfile = open(options['save'], 'wb')
    config.write(configfile)

# load results from file
# ------------------------------------------------------------------------------

def loadResult():
    if options['load']:
        config_parser = ConfigParser.RawConfigParser()
        config_parser.read(options['load'])
        if not config_parser.has_section('Sampling Result'):
            raise IOError("Invalid configuration file.")

        counts    = config.readMatrix(config_parser, 'Sampling Result', 'counts',    int)
        moments   = config.readMatrix(config_parser, 'Sampling Result', 'moments',   float)
        marginals = config.readMatrix(config_parser, 'Sampling Result', 'marginals', float)
        samples   = config.readVector(config_parser, 'Sampling Result', 'samples',   int)
        mpost     = config.readVector(config_parser, 'Sampling Result', 'mpost',     float)
        if config_parser.has_option('Sampling Result', 'bprob'):
            bprob = config.readVector(config_parser, 'Sampling Result', 'bprob',     float)
        else:
            bprob     = []
        result = {
            'moments'   : moments,
            'marginals' : marginals,
            'bprob'     : bprob,
            'mpost'     : mpost,
            'counts'    : counts,
            'samples'   : samples,
            'utility'   : [],
            'multibin_entropy'  : [],
            'differential_gain' : [] }
    else:
        result = {
            'moments'   : [],
            'marginals' : [],
            'bprob'     : [],
            'mpost'     : [],
            'counts'    : [],
            'samples'   : [],
            'utility'   : [],
            'multibin_entropy'  : [],
            'differential_gain' : [] }
    return result

# binning
# ------------------------------------------------------------------------------

def bin(counts, alpha, mprior):
    """Call the binning library."""
    counts_i = [ map(int, row) for row in counts ]
    alpha_i  = [ map(int, row) for row in alpha  ]
    mprior_i =   map(float, mprior)
    return interface.binning(counts_i, alpha_i, mprior_i, options)

# sampling
# ------------------------------------------------------------------------------

def selectItem(gain, counts, result):
    if options['strategy'] == 'uniform':
        return selectRandom(argmin(counts))
    elif options['strategy'] == 'uniform-random':
        return selectRandom(range(0,len(gain)))
    elif options['strategy'] == 'differential-gain':
        return selectRandom(argmax(gain))
    elif options['strategy'] == 'effective-counts':
        return selectRandom(argmin(result['effective_counts']))
    elif options['strategy'] == 'variance':
        stddev = map(math.sqrt, statistics.centralMoments(result['moments'], 2))
        return selectRandom(argmax(stddev))
    else:
        raise IOError('Unknown strategy: '+options['strategy'])

def experiment(ground_truth, index):
    if ground_truth[index] >= random.uniform(0.0, 1.0):
        return 0 # success
    else:
        return 1 # failure

def sampleFromGroundTruth(ground_truth, result, alpha, mprior):
    n = len(ground_truth)
    utility  = []
    marginal = options['marginal']
    options['marginal'] = 0
    if options['strategy'] == 'differential-gain':
        options['differential_gain'] = True
    if options['strategy'] == 'effective-counts':
        options['effective_counts'] = True
    if options['strategy'] == 'variance':
        options['n_moments'] = 3
    if result['counts']:
        counts = result['counts']
    else:
        counts = [ list(np.repeat(0, n)), list(np.repeat(0, n)) ]
    if result['samples']:
        samples = result['samples']
    else:
        samples = []
    for i in range(0, options['samples']):
        print "Sampling... %.1f%%" % ((float(i)+1)/float(options['samples'])*100)
        result  = bin(counts, alpha, mprior)
        gain    = map(lambda x: round(x, 4), result['differential_gain'])
        utility.append(gain[:])
        for j in range(0, options['blocks']):
            index   = selectItem(gain, map(sum, zip(*counts)), result)
            stddev = map(math.sqrt, statistics.centralMoments(result['moments'], 2))
            event   = experiment(ground_truth, index)
            samples.append(index)
            counts[event][index] += 1
            gain.pop(index)
    options['model_posterior'] = True
    options['n_moments'] = 3
    options['marginal'] = marginal
    result = bin(counts, alpha, mprior)
    result['counts']  = counts
    result['samples'] = samples
    result['utility'] = utility
    return result

# parse config
# ------------------------------------------------------------------------------

def parseConfig(config_file):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)

    if config_parser.sections() == []:
        raise IOError("Invalid configuration file.")
    if config_parser.has_section('Ground Truth'):
        gt     = config.readVector(config_parser, 'Ground Truth', 'gt', float)
        alpha  = config.readAlpha(config_parser, 2, len(gt), 'Ground Truth', int)
        mprior = config.readModelPrior(config_parser, len(gt), 'Ground Truth', int)
        options['script'] = config.readScript(config_parser, 'Ground Truth', os.path.dirname(config_file))
        result = loadResult()
        result = sampleFromGroundTruth(gt, result, alpha, mprior)
        if options['save']:
            saveResult(result)
        else:
            vis.plotSampling(result, gt, options)
            if options['plot-utility']:
                vis.plotUtilitySeries(result, options)
            show()

# main
# ------------------------------------------------------------------------------

options = {
    'blocks'            : 1,
    'samples'           : 0,
    'epsilon'           : 0.00001,
    'n_moments'         : 0,
    'marginal'          : 0,
    'marginal_step'     : 0.01,
    'marginal_range'    : (0.0,1.0),
    'which'             : 0,
    'strategy'          : 'differential-gain',
    'script'            : None,
    'load'              : None,
    'save'              : None,
    'plot-utility'      : False,
    'verbose'           : False,
    'prombsTest'        : False,
    'compare'           : False,
    'bprob'             : False,
    'differential_gain' : False,
    'effective_counts'  : False,
    'multibin_entropy'  : False,
    'model_posterior'   : False,
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=", "marginal", "marginal-range=",
                      "marginal-step=", "which=", "epsilon=", "moments", "blocks=",
                      "plot-utility", "strategy="]
        opts, tail = getopt.getopt(sys.argv[1:], "demr:s:k:n:bhvt", longopts)
    except getopt.GetoptError:
        usage()
        return 2
    output = None
    for o, a in opts:
        if o == "--blocks":
            options["blocks"] = int(a)
        if o in ("-v", "--verbose"):
            sys.stderr.write("Verbose mode turned on.\n")
            options["verbose"] = True
        if o in ("-t", "--prombsTest"):
            sys.stderr.write("Testing prombs.\n")
            options["prombsTest"] = True
        if o == "-n":
            options["samples"] = int(a)
        if o == "--strategy":
            options["strategy"] = a
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
        if o == "--plot-utility":
            options["plot-utility"] = True
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
