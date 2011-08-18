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
import operator
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
    import visualization as vis

import config
import interface
import statistics

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
    print "   -b                                 - compute break probabilities"
    print "   -d                                 - compute differential gain"
    print "       --lapsing=p                    - specify a lapsing probability"
    print "       --look-ahead=N                 - recursion depth for the sampling look ahead"
    print "   -m  --marginal                     - compute full marginal distribution"
    print "   -r  --marginal-range=(FROM,TO)     - limit range for the marginal distribution"
    print "   -s  --marginal-step=STEP           - step size for the marginal distribution"
    print "       --no-model-posterior           - do not compute the model posterior"
    print "       --epsilon=EPSILON              - epsilon for entropy estimations"
    print "   -n  --samples=N                    - number of samples"
    print "   -k  --moments=N                    - compute the first N>=2 moments"
    print "       --strategy=STRATEGY            - uniform, uniform-random, differential-gain (default),"
    print "                                        effective-counts, or variance"
    print "       --which=EVENT                  - for which event to compute the binning"
    print "       --algorithm=NAME               - select an algorithm [mgs, prombstree, default: prombs]"
    print "       --mgs-samples=BURN_IN:SAMPLES  - number of samples [default: 100:2000] for mgs"
    print
    print "       --port=PORT                    - connect to port from a matlab server for data collection"
    print
    print "       --threads=THREADS              - number of threads [default: 1]"
    print "       --stacksize=BYTES              - thread stack size [default: 256*1024]"
    print
    print "       --load=FILE                    - load result from file"
    print "       --save=FILE                    - save result to file"
    print "       --savefig=FILE                 - save figure to file"
    print
    print "   -h, --help                         - print help"
    print "   -v, --verbose                      - be verbose"
    print "   -t, --prombsTest                   - test prombs algorithm"
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
    configfile.close()

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
            'samples'   : samples }
    else:
        result = {
            'moments'   : [],
            'marginals' : [],
            'bprob'     : [],
            'mpost'     : [],
            'counts'    : [],
            'samples'   : [] }
    return result

# binning
# ------------------------------------------------------------------------------

def bin(counts_v, data, bin_options):
    """Call the binning library."""
    events = len(counts_v)
    counts = statistics.countStatistic(counts_v)
    alpha  = data['alpha']
    beta   = data['beta']
    gamma  = data['gamma']
    return interface.binning(events, counts, alpha, beta, gamma, bin_options)

# utility
# ------------------------------------------------------------------------------

hashutil = {}
hashexp  = {}

def prombsUtility(counts, data):
    bin_options = options.copy()
    bin_options['model_posterior'] = False
    bin_options['marginal']  = 0
    bin_options['n_moments'] = 0
    if options['strategy'] == 'variance':
        bin_options['n_moments'] = 2
    result = bin(counts, data, bin_options)
    if options['strategy'] == 'variance':
        return map(math.sqrt, statistics.centralMoments(result['moments'], 2))
    if options['strategy'] == 'differential-gain':
        return result['differential_gain']
    if options['strategy'] == 'effective-counts':
        return map(operator.neg, result['effective_counts'])

def prombsExpectation(y, counts, data):
    bin_options = options.copy()
    bin_options['differential_gain'] = False
    bin_options['effective_counts']  = False
    bin_options['model_posterior']   = False
    bin_options['marginal']  = 0
    bin_options['n_moments'] = 1
    bin_options['which'] = y
    result = bin(counts, data, bin_options)
    return result['moments'][0]

def computeKey(counts):
    return tuple(map(tuple, counts))

def recursiveUtility(counts, data, m, hashtree, hashutil, hashexp):
    key   = computeKey(counts)
    value = hashtree.get(key)
    if value:
        return value
    elif m == 0:
        utility = hashutil.get(key)
        if not utility:
            utility = prombsUtility(counts, data)
            hashutil[key] = utility
        return utility
    else:
        result = []
        if not hashexp.get(key):
            hashexp[key] = [ prombsExpectation(y, counts, data) for y in range(0, data['K']) ]
        for y in range(0, data['K']):
            utility     = []
            # ( pi(Y = y | X = y, X_n, Y_n) )_{x in X}
            expectation = hashexp.get(key)[y]
            for x in range(0, data['L']):
                counts[y][x] += 1
                tmp = recursiveUtility(counts, data, m-1, hashtree, hashutil, hashexp)
                counts[y][x] -= 1
                utility.append(max(tmp))
            result.append([ e*u for e, u in zip(expectation, utility) ])
        result = map(sum, zip(*result))
        hashtree[key] = result
        return result

def computeUtility(counts, data):
    return recursiveUtility(counts, data, options['look_ahead'], {}, hashutil, hashexp)

def selectItem(counts, data):
    if options['filter']:
        gainFilter = None
        exec options['filter']
        if not gainFilter is None:
            # TODO
            gain = gainFilter(gain, map(sum, zip(*counts)), result)
    if options['strategy'] == 'uniform':
        utility = map(operator.neg, map(sum, zip(*counts)))
    elif options['strategy'] == 'uniform-random':
        utility = [ 0.0 for i in range(0,len(gain)) ]
    elif options['strategy'] == 'differential-gain':
        utility = computeUtility(counts, data)
    elif options['strategy'] == 'effective-counts':
        utility = computeUtility(counts, data)
    elif options['strategy'] == 'variance':
        utility = computeUtility(counts, data)
    else:
        raise IOError('Unknown strategy: '+options['strategy'])
    return selectRandom(argmax(utility)), utility

# experiment
# ------------------------------------------------------------------------------

def experiment(index, data, msocket):
    if data['gt']:
        # load previous state of the random generator
        if data['rand_states']:
            random.setstate(data['rand_states'][index])
        # lapsing
        if options['lapsing'] >= random.uniform(0.0, 1.0):
            if 0.5 >= random.uniform(0.0, 1.0):
                return 0 # success
            else:
                return 1 # failure
        # real psychometric function
        sample = random.uniform(0.0, 1.0)
        # save state of the random generator
        if data['rand_states']:
            data['rand_states'][index] = random.getstate()
        if data['gt'][index] >= sample:
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

# sampling
# ------------------------------------------------------------------------------

def sample(result, data):
    msocket   = None
    if options['port']:
        msocket = open_msocket()
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
        index, utility = selectItem(counts, data)
        event  = experiment(index, data, msocket)
        samples.append(index)
        counts[event][index] += 1
    index, utility = selectItem(counts, data)
    result = bin(counts, data, options)
    result['counts']  = counts
    result['samples'] = samples
    result['differential_gain'] = utility
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
        config.readVisualization(config_parser, 'Ground Truth', os.path.dirname(config_file), options)
        config.readFilter(config_parser, 'Ground Truth', os.path.dirname(config_file), options)
        config.readAlgorithm(config_parser, 'Ground Truth', os.path.dirname(config_file), options)
        config.readMgsSamples(config_parser, 'Ground Truth', os.path.dirname(config_file), options)
        data['rand_states'] = config.readSeeds(config_parser, 'Ground Truth', 'seeds')
        data['gt'] = config.readVector(config_parser, 'Ground Truth', 'gt', float)
        data['L'] = len(data['gt'])
        data['alpha'], data['beta'], data['gamma'] = \
            config.getParameters(config_parser, 'Ground Truth', os.path.dirname(config_file), data['K'], data['L'])
        config.readStrategy(config_parser, 'Ground Truth', options)
        result = loadResult()
        result = sample(result, data)
    if config_parser.has_section('Experiment'):
        config.readVisualization(config_parser, 'Experiment', os.path.dirname(config_file), options)
        config.readFilter(config_parser, 'Experiment', os.path.dirname(config_file), options)
        config.readAlgorithm(config_parser, 'Experiment', os.path.dirname(config_file), options)
        config.readMgsSamples(config_parser, 'Experiment', os.path.dirname(config_file), options)
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
            show()

# main
# ------------------------------------------------------------------------------

options = {
    'samples'           : 0,
    'look_ahead'        : 0,
    'epsilon'           : 0.00001,
    'n_moments'         : 2,
    'mgs_samples'       : (100,2000),
    'marginal'          : 0,
    'marginal_step'     : 0.01,
    'marginal_range'    : (0.0,1.0),
    'which'             : 0,
    'lapsing'           : 0.0,
    'threads'           : 1,
    'stacksize'         : 256*1024,
    'algorithm'         : 'prombs',
    'strategy'          : 'differential-gain',
    'port'              : None,
    'filter'            : None,
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
        longopts   = ["help", "verbose", "load=", "save=", "marginal", "marginal-range=",
                      "marginal-step=", "which=", "epsilon=", "moments", "look-ahead=",
                      "strategy=", "savefig=", "lapsing=", "port=", "threads=", "stacksize=",
                      "algorithm=", "samples=", "mgs-samples", "no-model-posterior"]
        opts, tail = getopt.getopt(sys.argv[1:], "dmr:s:k:n:bhvt", longopts)
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
        if o == "--look-ahead":
            options["look_ahead"] = int(a)
        if o in ("-n", "--samples"):
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
            if int(a) >= 2:
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
        if o == "--algorithm":
            options["algorithm"] = a
        if o == "--mgs-samples":
            options["mgs_samples"] = tuple(map(int, a.split(":")))
        if o == "--no-model-posterior":
            options["model_posterior"] = False
    if options['strategy'] == 'differential-gain':
        options['differential_gain'] = True
    if options['strategy'] == 'effective-counts':
        options['effective_counts'] = True
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
