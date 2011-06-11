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
import socket
import random

from itertools import izip

def importMatplotlib(backend=None):
    global vis
    from matplotlib import use
    if backend:
        use(backend)
    import bayesian_binning.visualization as vis

import bayesian_binning.config        as config
import bayesian_binning.interface     as interface
import bayesian_binning.statistics    as statistics

# global options
# ------------------------------------------------------------------------------

verbose = False

# usage
# ------------------------------------------------------------------------------

def usage():
    """Print usage."""
    print
    print "adaptive-sampling [option]... FILE "
    print
    print "Options:"
    print "   -b                             - compute break probabilities"
    print "       --blocks=N                 - sample in blocks of N measurements"
    print "   -d                             - compute differential gain"
    print "       --lapsing=p                - specify a lapsing probability"
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
    print "       --port=PORT                - connect to port from a matlab server for data collection"
    print
    print "       --threads=THREADS          - number of threads [default: 1]"
    print "       --stacksize=BYTES          - thread stack size [default: 256*1024]"
    print
    print "       --load=FILE                - load result from file"
    print "       --save=FILE                - save result to file"
    print "       --savefig=FILE             - save figure to file"
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

def close_msocket(m_socket):
    m_socket.send('9999\n')
    m_socket.shutdown(1)
    m_socket.close()

def open_msocket():
    msocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    msocket.connect(('localhost', options['port']))
    if verbose:
        print "Connected to localhost on port: " + str(options['port'])
    return msocket

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
            'differential_gain' : [] }
    return result

# binning
# ------------------------------------------------------------------------------

def bin(counts_v, data):
    """Call the binning library."""
    events = len(counts_v)
    counts = statistics.countStatistic(counts_v)
    alpha  = data['alpha']
    beta   = data['beta']
    gamma  = data['gamma']
    return interface.binning(events, counts, alpha, beta, gamma, options)

# sampling
# ------------------------------------------------------------------------------

def selectItem(gain, counts, result):
    if options['filter']:
        gainFilter = None
        exec options['filter']
        if not gainFilter is None:
            gain = gainFilter(gain, counts, result)
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

def experiment(index, data, msocket):
    if data['gt']:
        # lapsing
        if options['lapsing'] >= random.uniform(0.0, 1.0):
            if 0.5 >= random.uniform(0.0, 1.0):
                return 0 # success
            else:
                return 1 # failure
        # real psychometric function
        if data['gt'][index] >= random.uniform(0.0, 1.0):
            return 0 # success
        else:
            return 1 # failure
    else:
        if options['port']:
            msocket.send(str(index) + '\n')
            ret = None
            while not ret == 0 and not ret == 1:
                ret = int(msocket.recv(1))
            if verbose:
                print str(index) + ' answer: ' + str(ret)
        else:
            print index
            ret = None
            while not ret == 0 and not ret == 1:
                ret = int(raw_input())
    return ret

def sample(result, data):
    utility  = []
    marginal = options['marginal']
    msocket  = None
    options['marginal'] = 0
    if options['port']:
        msocket = open_msocket()
    if options['strategy'] == 'differential-gain':
        options['differential_gain'] = True
    if options['strategy'] == 'effective-counts':
        options['effective_counts'] = True
    if options['strategy'] == 'variance':
        options['n_moments'] = 3
    if result['counts']:
        counts = result['counts']
    else:
        counts = [ list(np.repeat(0, data['L'])),
                   list(np.repeat(0, data['L'])) ]
    if result['samples']:
        samples = result['samples']
    else:
        samples = []
    for i in range(0, options['samples']):
        print >> sys.stderr, "Sampling... %.1f%%" % ((float(i)+1)/float(options['samples'])*100)
        result  = bin(counts, data)
        gain    = map(lambda x: round(x, 4), result['differential_gain'])
        utility.append(gain[:])
        for j in range(0, options['blocks']):
            index  = selectItem(gain, map(sum, zip(*counts)), result)
            stddev = map(math.sqrt, statistics.centralMoments(result['moments'], 2))
            event  = experiment(index, data, msocket)
            samples.append(index)
            counts[event][index] += 1
            gain.pop(index)
    options['model_posterior'] = True
    options['n_moments'] = 3
    options['marginal'] = marginal
    result = bin(counts, data)
    result['counts']  = counts
    result['samples'] = samples
    result['utility'] = utility
    if options['port']:
        close_msocket(msocket)
    return result

# parse config
# ------------------------------------------------------------------------------

def parseConfig(config_file):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)

    data = { 'K'     : 2,
             'L'     : 0,
             'gt'    : None,
             'alpha' : None,
             'beta'  : None,
             'gamma' : None }

    if config_parser.sections() == []:
        raise IOError("Invalid configuration file.")
    if config_parser.has_section('Ground Truth'):
        options['visualization'] = config.readVisualization(config_parser, 'Ground Truth', os.path.dirname(config_file))
        options['filter'] = config.readFilter(config_parser, 'Ground Truth', os.path.dirname(config_file))
        data['gt']  = config.readVector(config_parser, 'Ground Truth', 'gt', float)
        data['L']   = len(data['gt'])
        data['alpha'], data['beta'], data['gamma'] = \
            config.getParameters(config_parser, 'Ground Truth', os.path.dirname(config_file), data['K'], data['L'])
        config.readStrategy(config_parser, 'Ground Truth', options)
        result = loadResult()
        result = sample(result, data)
    if config_parser.has_section('Experiment'):
        options['visualization'] = config.readVisualization(config_parser, 'Experiment', os.path.dirname(config_file))
        options['filter'] = config.readFilter(config_parser, 'Experiment', os.path.dirname(config_file))
        data['L'] = int(config_parser.get('Experiment', 'bins'))
        data['alpha'], data['beta'], data['gamma'] = \
            config.getParameters(config_parser, 'Experiment', os.path.dirname(config_file), data['K'], data['L'])
        config.readStrategy(config_parser, 'Experiment', options)
        result = loadResult()
        result = sample(result, data)
    if options['save']:
        saveResult(result)
    else:
        if options['savefig']:
            importMatplotlib('Agg')
            from matplotlib.pyplot import savefig
            vis.plotSampling(result, options, data)
            savefig(options['savefig'], bbox_inches='tight', pad_inches=0)
        else:
            importMatplotlib()
            from matplotlib.pyplot import show
            vis.plotSampling(result, options, data)
            if options['plot-utility']:
                vis.plotUtilitySeries(result, options, data)
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
    'lapsing'           : 0.0,
    'threads'           : 1,
    'stacksize'         : 256*1024,
    'strategy'          : 'differential-gain',
    'port'              : None,
    'filter'            : None,
    'visualization'     : None,
    'load'              : None,
    'save'              : None,
    'savefig'           : None,
    'plot-utility'      : False,
    'verbose'           : False,
    'prombsTest'        : False,
    'compare'           : False,
    'bprob'             : False,
    'differential_gain' : False,
    'effective_counts'  : False,
    'model_posterior'   : False,
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=", "marginal", "marginal-range=",
                      "marginal-step=", "which=", "epsilon=", "moments", "blocks=",
                      "plot-utility", "strategy=", "savefig=", "lapsing=", "port=",
                      "threads=", "stacksize=" ]
        opts, tail = getopt.getopt(sys.argv[1:], "dmr:s:k:n:bhvt", longopts)
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
        if o == "--lapsing":
            options["lapsing"] = float(a)
        if o == "-d":
            options["differential_gain"] = True
        if o == "--load":
            options["load"] = a
        if o == "--port":
            options["port"] = int(a)
        if o == "--save":
            options["save"] = a
        if o == "--savefig":
            options["savefig"] = a
        if o == "--which":
            options["which"] = int(a)
        if o == "--epsilon":
            options["epsilon"] = float(a)
        if o == "--plot-utility":
            options["plot-utility"] = True
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
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
