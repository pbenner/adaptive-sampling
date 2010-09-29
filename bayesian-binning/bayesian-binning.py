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

import getopt
import os
import interface
import ConfigParser
import numpy as np
import math
from matplotlib import *
from matplotlib.pyplot import *

# global options
# ------------------------------------------------------------------------------

verbose = False

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print "bayesian-binning.py [option]... FILE "
    print
    print "Options:"
    print "   -h, --help         - print help"
    print "   -v, --verbose      - be verbose"

# plotting
# ------------------------------------------------------------------------------

def plotspikes(ax, x, timings):
    """Plot trials of spike trains."""
    X = []
    Y = []
    for i in range(0, len(timings)):
        for val in timings[i]:
            X.append(val)
            Y.append(i)
    ax.set_xlim(x[0],x[-1])
    ax.set_ylabel('Trials')
    ax.plot(X, Y, 'k|')

def plotbin(ax, x, result):
    """Plot the binning result."""
    font = {'family'     : 'serif',
            'color'      : 'k',
            'weight'     : 'normal',
            'size'       : 12 }

    N = len(result[0])
    if x==None:
        x = np.arange(0, N, 1)
    ax.set_xlim(x[0],x[-1])
    ax.plot(x, [a + b for a, b in zip(result[0], result[1])], 'k--')
    ax.plot(x, [a - b for a, b in zip(result[0], result[1])], 'k--')
    ax.plot(x, result[0], 'r')
    ax.set_xlabel('bin', font)
    ax.set_ylabel('P(x)', font)

# binning
# ------------------------------------------------------------------------------

def bin(counts, trials, prior):
    """Call the binning library."""
    return interface.binning(counts, trials, prior)

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

def parseConfig(file):
    config = ConfigParser.RawConfigParser()
    config.read(file)

    if config.has_section('Counts'):
        trials      = config.getint('Counts', 'trials')
        counts_str  = config.get   ('Counts', 'counts')
        counts      = []

        for value in counts_str.split(' '):
            counts.append(int(value))

        N           = len(counts)
        prior       = list(np.repeat(1, N))
        #prior       = list(np.repeat(0, N))
        #prior[1]    = 1
        fig = figure()
        ax1 = fig.add_subplot(1,2,1)
        plotbin(ax1, None, bin(counts, trials, prior))
        show()
    if config.has_section('Trials'):
        binsize     = config.getint('Trials', 'binsize')
        timings_str = config.get   ('Trials', 'timings')
        timings     = []
        for line in timings_str.split('\n'):
            if line != '':
                timings.append([int(a) for a in line.split(' ')])

        x, counts   = timingsToCounts(timings, binsize)
        trials      = len(timings)
        N           = len(counts)
        prior       = list(np.repeat(1, N))
        fig = figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        plotspikes(ax1, x, timings)
        plotbin(ax2, x, bin(counts, trials, prior))
        show()

# main
# ------------------------------------------------------------------------------

def main():
    global verbose
    try:
        opts, tail = getopt.getopt(sys.argv[1:], "hv", ["help", "output="])
    except getopt.GetoptError:
        usage()
        return 2
    output = None
    for o, a in opts:
        if o == "-v":
            print "verbose mode turned on"
            verbose = True
        if o in ("-h", "--help"):
            usage()
            return 0
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
