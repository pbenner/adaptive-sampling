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
import sys
import ConfigParser
import numpy as np
import math

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "entropy.py [option]... FILE "
    print
    print "Options:"
    print
    print "   -h, --help                  - print help"
    print

# usage
# ------------------------------------------------------------------------------

def computeEntropy(file):
    config = ConfigParser.RawConfigParser()
    config.read(file)
    if not config.has_section('Result'):
        raise IOError("Invalid configuration file.")
    pdf_str     = config.get('Result', 'pdf')
    pdf         = map(float, pdf_str.split(' '))
    entropy     = 0
    for p in pdf:
        entropy += p*math.log(p)
        entropy += (1-p)*math.log(1-p)
    print -entropy

def main():
    global options
    try:
        longopts   = ["help"]
        opts, tail = getopt.getopt(sys.argv[1:], "h", longopts)
    except getopt.GetoptError:
        usage()
        return 2
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            return 0
    if len(tail) != 1:
        usage()
        return 1
    computeEntropy(tail[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
