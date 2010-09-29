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

# binning
# ------------------------------------------------------------------------------

def bin(file):
    """Read config file and call the binning library."""
    config = ConfigParser.RawConfigParser()
    config.read(file)

    trials     = config.getint('Data', 'trials')
    counts_str = config.get   ('Data', 'counts')
    counts     = []

    for value in counts_str.split(' '):
        counts.append(int(value))

    N          = len(counts)
    prior      = list(np.repeat(1, N))
    #prior      = list(np.repeat(0, N))
    #prior[1]   = 1
    result     = interface.binning(counts, trials, prior)
    return result

# plotting
# ------------------------------------------------------------------------------

def plotbin(result):
    """Plot the binning result."""
    font = {'family'     : 'serif',
            'color'      : 'k',
            'weight'     : 'normal',
            'size'       : 12 }

    N = len(result[0])
    x = np.arange(0, N, 1)
    plot(x, [a + b for a, b in zip(result[0], result[1])], 'k--')
    plot(x, [a - b for a, b in zip(result[0], result[1])], 'k--')
    plot(x, result[0], 'r')
    xlabel('bin', font)
    ylabel('P(x)', font)
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
    result = bin(tail[0])
    plotbin(result)
    return 0

if __name__ == "__main__":
    sys.exit(main())
