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
import interface
import ConfigParser
import numpy as np
import math
from itertools import izip
from matplotlib import *
from matplotlib.pyplot import *
import matplotlib.patches as patches
import matplotlib.path as path

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
    print "       --likelihood=LIKELIHOOD - multinomial, binomial"
    print "   -s, --sigma=SIGMA           - sigma hyperparameter"
    print "   -g, --gamma=GAMMA           - gamma hpyerparameter"
    print "   -b                          - compute break probabilities"
    print
    print "   -m                          - use GNU multiple precision library"
    print "       --compare               - compare log scale with gmp library"
    print
    print "       --load                  - load result from file"
    print "       --save                  - save result to file"
    print
    print "   -h, --help                  - print help"
    print "   -v, --verbose               - be verbose"
    print

# tools
# ------------------------------------------------------------------------------

argmax = lambda array:max(izip(array, xrange(len(array))))

# plotting
# ------------------------------------------------------------------------------
font = {'family'     : 'serif',
        'color'      : 'k',
        'weight'     : 'normal',
        'size'       : 12 }

def plotmodelpost(ax, modelpost):
    N = len(modelpost)
    x = np.arange(0, N+1, 1)
    modelpost.insert(0, 0)
    ax.step(x, modelpost, 'r--', where='mid', linewidth=1)
    ax.grid(True)

    left    = np.array(x[:-1]) + 0.5
    right   = np.array(x[1:])  + 0.5
    bottom  = np.zeros(len(left))
    top     = bottom + modelpost[1:]
    XY      = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    barpath = path.Path.make_compound_path_from_polys(XY)
    patch   = patches.PathPatch(barpath, facecolor='green', edgecolor='gray', alpha=0.8)
    ax.add_patch(patch)
    ax.set_xlabel(r'$M$',  font)
    ax.set_ylabel(r'$P(M|D)$', font)

def plotspikes(ax, x, timings):
    """Plot trials of spike trains."""
    X = []
    Y = []
    for i in range(0, len(timings)):
        for val in timings[i]:
            X.append(val)
            Y.append(i)
    ax.set_xlim(x[0],x[-1])
    ax.set_ylabel(r'$t$')
    ax.set_ylabel('Trial')
    ax.plot(X, Y, 'k|')

def plotbinboundaries(ax, x, bprob, modelpost):
    ax.plot(x[1:-1], bprob[1:-1], 'g')
    ax.set_ylim(0,1)
    ax.set_ylabel('P(Break|D)', font)
#    nbins = argmax(modelpost)[1]+1
#    bprob_max = sorted(bprob)[-nbins:]
#    for b in bprob_max:
#        i = bprob.index(b)
#        ax.axvline(x[i])

def plotbin(ax, x, exp, var, bprob, modelpost):
    """Plot the binning result."""
    N = len(exp)
    if x==None:
        x = np.arange(0, N, 1)
    ax.set_xlim(x[0],x[-1])
    ax.plot(x, [min(1,a + math.sqrt(b)) for a, b in zip(exp, var)], 'k--')
    ax.plot(x, [max(0,a - math.sqrt(b)) for a, b in zip(exp, var)], 'k--')
    ax.plot(x, exp, 'r')
    ax.set_xlabel('t',  font)
    ax.set_ylabel(r'$P(S_i|D)$', font)
    if bprob:
        plotbinboundaries(ax.twinx(), x, bprob, modelpost)

# binning
# ------------------------------------------------------------------------------

def bin(successes, failures, mprior):
    """Call the binning library."""
    if options['load']:
        config = ConfigParser.RawConfigParser()
        config.read(options['load'])
        if not config.has_section('Result'):
            raise IOError("Invalid configuration file.")
        pdf_str   = config.get('Result', 'pdf')
        var_str   = config.get('Result', 'var')
        mpost_str = config.get('Result', 'mpost')
        pdf       = map(float, pdf_str.split(' '))
        var       = map(float, var_str.split(' '))
        mpost     = map(float, mpost_str.split(' '))
        if config.has_option('Result', 'bprob'):
            bprob_str = config.get('Result', 'bprob')
            bprob     = map(float, bprob_str.split(' '))
        else:
            bprob = []
        return [pdf, var, bprob, mpost]
    elif options["compare"]:
        successes_i = map(int, successes)
        failures_i  = map(int, failures)
        mprior_i    = map(float, mprior)
        options["gmp"] = False
        result1 = interface.binning(successes_i, failures_i, mprior_i, options)
        options["gmp"] = True
        result2 = interface.binning(successes_i, failures_i, mprior_i, options)
        print [abs(a - b) for a, b in zip(result1[0], result2[0])]
        exit(0)
    else:
        successes_i = map(int, successes)
        failures_i  = map(int, failures)
        mprior_i    = map(float, mprior)
        alpha_i     = map(int, [options['sigma'], options['gamma']])
        return interface.binning([successes_i, failures_i], alpha_i, mprior_i, options)

# save result
# ------------------------------------------------------------------------------

def saveResult(result):
    config = ConfigParser.ConfigParser()
    config.add_section('Result')
    config.set('Result', 'pdf', " ".join(map(str, result[0])))
    config.set('Result', 'var', " ".join(map(str, result[1])))
    config.set('Result', 'bprob', " ".join(map(str, result[2])))
    config.set('Result', 'mpost', " ".join(map(str, result[3])))
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
        raise ValueError("number of trials is smaller than some counts")
    return failures

def readMPrior(models_str, N):
    models = []
    mprior = list(np.repeat(0, N))
    for str in models_str.split(' '):
        models.append(int(str))
    num_models = len(models)
    for model in models:
        mprior[model-1] = 1.0/num_models
    return mprior

def readOptions(config, section):
    if config.has_option(section, 'likelihood'):
        likelihood_str = config.get(section, 'likelihood')
        if likelihood_str == "multinomial":
            options["likelihood"] = 1
        if likelihood_str == "binomial":
            options["likelihood"] = 2
    if config.has_option(section, 'sigma'):
        options["sigma"] = int(config.get(section, 'sigma'))
    if config.has_option(section, 'gamma'):
        options["gamma"] = int(config.get(section, 'gamma'))

def parseConfig(file):
    config = ConfigParser.RawConfigParser()
    config.read(file)

    if config.sections() == []:
        raise IOError("Invalid configuration file.")
    if config.has_section('Counts'):
        readOptions(config, 'Counts')
        trials      = config.getint('Counts', 'trials')
        counts_str  = config.get   ('Counts', 'counts')
        successes   = []
        for value in counts_str.split(' '):
            successes.append(int(value))
        failures    = computeFailures(successes, trials)

        N           = len(successes)
        prior       = list(np.repeat(1, N))
        if config.has_option('Counts', 'mprior'):
            prior   = readMPrior(config.get('Counts', 'mprior'), N)
        result      = bin(successes, failures, prior)
        if options['save']:
            saveResult(result)
        else:
            fig = figure()
            fig.subplots_adjust(hspace=0.35)
            ax1 = fig.add_subplot(2,1,1)
            ax2 = fig.add_subplot(2,1,2)
            plotbin(ax1, None, result[0], result[1], result[2], result[3])
            plotmodelpost(ax2, result[3])
            show()
    if config.has_section('Trials'):
        readOptions(config, 'Trials')
        binsize     = config.getint('Trials', 'binsize')
        timings_str = config.get   ('Trials', 'timings')
        timings     = []
        for line in timings_str.split('\n'):
            if line != '':
                timings.append([int(a) for a in line.split(' ')])

        x, successes = timingsToCounts(timings, binsize)
        trials       = len(timings)
        failures     = computeFailures(successes, trials)
        N            = len(successes)
        prior        = list(np.repeat(1, N))
        if config.has_option('Trials', 'mprior'):
            prior    = readMPrior(config.get('Trials', 'mprior'), N)
        result       = bin(successes, failures, prior)
        if options['save']:
            saveResult(result)
        else:
            fig = figure()
            fig.subplots_adjust(hspace=0.35)
            ax1 = fig.add_subplot(3,1,1)
            ax2 = fig.add_subplot(3,1,2)
            ax3 = fig.add_subplot(3,1,3)
            plotspikes(ax1, x, timings)
            plotbin   (ax2, x, result[0], result[1], result[2], result[3])
            plotmodelpost(ax3, result[3])
            show()

# main
# ------------------------------------------------------------------------------

options = {
    'verbose'    : False,
    'gmp'        : False,
    'compare'    : False,
    'bprob'      : False,
    'likelihood' : 1,
    'sigma'      : 1,
    'gamma'      : 1,
    'load'       : None,
    'save'       : None
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "compare", "likelihood=",
                      "sigma=", "gamma=", "load=", "save="]
        opts, tail = getopt.getopt(sys.argv[1:], "bhvs:g:m", longopts)
    except getopt.GetoptError:
        usage()
        return 2
    output = None
    for o, a in opts:
        if o in ("-v", "--verbose"):
            sys.stderr.write("Verbose mode turned on.\n")
            options["verbose"] = True
        if o in ("-h", "--help"):
            usage()
            return 0
        if o == "-b":
            options["bprob"] = True
        if o == "-m":
            sys.stderr.write("Using GNU multiple precision library.\n")
            options["gmp"] = True
        if o == "--compare":
            sys.stderr.write("Comparing log scale with GNU multiple precision library.\n")
            options["compare"] = True
        if o == "--likelihood":
            if a == "multinomial":
                options["likelihood"] = 1
            if a == "binomial":
                options["likelihood"] = 2
        if o in ("-s", "--sigma"):
            options["sigma"] = int(a)
        if o in ("-g", "--gamma"):
            options["gamma"] = int(a)
        if o == "--load":
            options["load"] = a
        if o == "--save":
            options["save"] = a
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
