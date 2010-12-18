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
    print "   -b                          - compute break probabilities"
    print "       --likelihood=LIKELIHOOD - multinomial, binomial"
    print "       --which=EVENT           - for which event to compute the binning"
    print "       --epsilon=EPSILON       - epsilon for entropy estimations"
    print
    print "       --load                  - load result from file"
    print "       --save                  - save result to file"
    print
    print "   -h, --help                  - print help"
    print "   -v, --verbose               - be verbose"
    print "   -t, --prombsTest            - test prombs algorithm"
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
    ax.set_xlabel(r'$m_B$',  font)
    ax.set_ylabel(r'$P(m_B|D)$', font)

def plotentropy(ax, entropy):
    N = len(entropy)
    x = np.arange(0, N+1, 1)
    entropy.insert(0, 0)
    ax.step(x, entropy, 'r--', where='mid', linewidth=1)
    ax.grid(True)

    left    = np.array(x[:-1]) + 0.5
    right   = np.array(x[1:])  + 0.5
    bottom  = np.zeros(len(left))
    top     = bottom + entropy[1:]
    XY      = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    barpath = path.Path.make_compound_path_from_polys(XY)
    patch   = patches.PathPatch(barpath, facecolor='green', edgecolor='gray', alpha=0.8)
    ax.add_patch(patch)
    ax.set_xlabel(r'$m_B$',  font)
    ax.set_ylabel(r'$H(\mathcal{P}(\mathcal{X})|D,m_B)$', font)

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
    if bprob and options['bprob']:
        plotbinboundaries(ax.twinx(), x, bprob, modelpost)

# binning
# ------------------------------------------------------------------------------

def bin(counts, alpha, mprior):
    """Call the binning library."""
    if options['load']:
        config = ConfigParser.RawConfigParser()
        config.read(options['load'])
        if not config.has_section('Result'):
            raise IOError("Invalid configuration file.")
        pdf_str     = config.get('Result', 'pdf')
        var_str     = config.get('Result', 'var')
        mpost_str   = config.get('Result', 'mpost')
        entropy_str = config.get('Result', 'entropy')
        pdf         = map(float, pdf_str.split(' '))
        var         = map(float, var_str.split(' '))
        mpost       = map(float, mpost_str.split(' '))
        entropy     = map(float, entropy_str.split(' '))
        if config.has_option('Result', 'bprob'):
            bprob_str = config.get('Result', 'bprob')
            bprob     = map(float, bprob_str.split(' '))
        else:
            bprob = []
        return [pdf, var, bprob, mpost, entropy]
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
    config.set('Result', 'pdf',     " ".join(map(str, result[0])))
    config.set('Result', 'var',     " ".join(map(str, result[1])))
    config.set('Result', 'bprob',   " ".join(map(str, result[2])))
    config.set('Result', 'mpost',   " ".join(map(str, result[3])))
    config.set('Result', 'entropy', " ".join(map(str, result[4])))
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
        readOptions(config, 'Counts')
        counts_str  = config.get   ('Counts', 'counts')
        counts      = []
        for line in counts_str.split('\n'):
            if line != '':
                counts.append([int(a) for a in line.split(' ')])
        N           = len(counts[0])
        prior       = list(np.repeat(1, N))
        alpha       = readAlpha(config, 'Counts', len(counts))
        if config.has_option('Counts', 'mprior'):
            prior   = readMPrior(config.get('Counts', 'mprior'), N)
        result      = bin(counts, alpha, prior)
        if options['save']:
            saveResult(result)
        else:
            fig = figure()
            fig.subplots_adjust(hspace=0.35)
            ax1 = fig.add_subplot(3,1,1)
            ax2 = fig.add_subplot(3,1,2)
            ax3 = fig.add_subplot(3,1,3)
            plotbin(ax1, None, result[0], result[1], result[2], result[3])
            plotmodelpost(ax2, result[3])
            plotentropy(ax3, result[4])
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
        counts       = [successes, failures]
        N            = len(successes)
        prior        = list(np.repeat(1, N))
        alpha        = readAlpha(config, 'Trials', 2)
        if config.has_option('Trials', 'mprior'):
            prior    = readMPrior(config.get('Trials', 'mprior'), N)
        result       = bin(counts, alpha, prior)
        if options['save']:
            saveResult(result)
        else:
            fig = figure()
            fig.subplots_adjust(hspace=0.35)
            ax1 = fig.add_subplot(4,1,1)
            ax2 = fig.add_subplot(4,1,2)
            ax3 = fig.add_subplot(4,1,3)
            ax4 = fig.add_subplot(4,1,4)
            plotspikes(ax1, x, timings)
            plotbin   (ax2, x, result[0], result[1], result[2], result[3])
            plotmodelpost(ax3, result[3])
            plotentropy(ax4, result[4])
            show()

# main
# ------------------------------------------------------------------------------

options = {
    'epsilon'    : 0.00001,
    'verbose'    : False,
    'prombsTest' : False,
    'compare'    : False,
    'bprob'      : False,
    'likelihood' : 1,
    'load'       : None,
    'save'       : None,
    'which'      : 0
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "likelihood=", "load=",
                      "save=", "which=", "epsilon=", "prombsTest"]
        opts, tail = getopt.getopt(sys.argv[1:], "bhvt", longopts)
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
        if o in ("-h", "--help"):
            usage()
            return 0
        if o == "-b":
            options["bprob"] = True
        if o == "--likelihood":
            if a == "multinomial":
                options["likelihood"] = 1
            if a == "binomial":
                options["likelihood"] = 2
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
