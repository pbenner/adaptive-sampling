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
import Queue
import threading

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
    print "       --lapsing=p                    - specify a lapsing probability"
    print "       --look-ahead=N                 - recursion depth for the sampling look ahead"
    print "   -m  --marginal                     - compute full marginal distribution"
    print "   -r  --marginal-range=(FROM,TO)     - limit range for the marginal distribution"
    print "   -s  --marginal-step=STEP           - step size for the marginal distribution"
    print "       --no-model-posterior           - do not compute the model posterior"
    print "       --epsilon=EPSILON              - epsilon for entropy estimations"
    print "   -n  --samples=N                    - number of samples"
    print "   -k  --moments=N                    - compute the first N>=2 moments"
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
    print "Sampling strategies:"
    print "       --strategy=STRATEGY            - uniform, uniform-random, entropy (default),"
    print "                                        effective-counts, or variance"
    print "       --differential-entropy         - include differential entropy for entropy based sampling"
    print "       --multibin-entropy             - include multibin entropy for entropy based sampling"
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
    config.set('Sampling Result', 'entropy',   " ".join(map(str, result['entropy'])))
    config.set('Sampling Result', 'bprob',     " ".join(map(str, result['bprob'])))
    config.set('Sampling Result', 'mpost',     " ".join(map(str, result['mpost'])))
    config.set('Sampling Result', 'states',    "\n".join(map(str, result['states'])))
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
        entropy   = config.readVector(config_parser, 'Sampling Result', 'entropy',   float)
        states    = config.readStates(config_parser, 'Sampling Result', 'states')
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
            'entropy'   : entropy,
            'states'    : states }
    else:
        result = {
            'moments'   : [],
            'marginals' : [],
            'bprob'     : [],
            'mpost'     : [],
            'counts'    : [],
            'samples'   : [],
            'entropy'   : [],
            'states'    : [] }
    return result

# binning interface
# ------------------------------------------------------------------------------

def bin_entropy(counts_v, data, bin_options):
    """Call the binning library."""
    events = len(counts_v)
    counts = statistics.countStatistic(counts_v)
    alpha  = data['alpha']
    beta   = data['beta']
    gamma  = data['gamma']
    return interface.entropy(events, counts, alpha, beta, gamma, bin_options)

def bin(counts_v, data, bin_options):
    """Call the binning library."""
    events = len(counts_v)
    counts = statistics.countStatistic(counts_v)
    alpha  = data['alpha']
    beta   = data['beta']
    gamma  = data['gamma']
    return interface.binning(events, counts, alpha, beta, gamma, bin_options)

# prombs wrapper
# ------------------------------------------------------------------------------

def prombsEntropy(counts, data):
    bin_options = options.copy()
    return bin_entropy(counts, data, bin_options)

def prombsUtility(counts, data):
    bin_options = options.copy()
    bin_options['utility']         = True
    bin_options['model_posterior'] = False
    bin_options['marginal']        = 0
    bin_options['n_moments']       = 0
    if options['strategy'] == 'variance':
        bin_options['n_moments'] = 2
    result = bin(counts, data, bin_options)
    if options['strategy'] == 'variance':
        return map(math.sqrt, statistics.centralMoments(result['moments'], 2))
    if options['strategy'] == 'entropy':
        return result['utility']
    if options['strategy'] == 'effective-counts':
        return result['utility']

def prombsExpectation(y, counts, data):
    bin_options = options.copy()
    bin_options['utility']         = False
    bin_options['model_posterior'] = False
    bin_options['marginal']        = 0
    bin_options['n_moments']       = 1
    bin_options['which']           = y
    result = bin(counts, data, bin_options)
    return result['moments'][0]

# recursive computation of utilities
# ------------------------------------------------------------------------------

hashutil = {}
hashexp  = {}

def computeKey(counts):
    return tuple(map(tuple, counts))

def recursiveUtility(counts, data, m, hashtree, hashutil, hashexp):
    key   = computeKey(counts)
    value = hashtree.get(key)
    if value:
        return value
    elif m == 0:
        if not hashutil.get(key):
            hashutil[key] = prombsUtility(counts, data)
        return hashutil[key]
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

# precompute utility
# ------------------------------------------------------------------------------

class ThreadUtility(threading.Thread):
    def __init__(self, queue_in, queue_out, data):
        threading.Thread.__init__(self)
        self.queue_in  = queue_in
        self.queue_out = queue_out
        self.data      = data
    def run(self):
        while True:
            counts = self.queue_in.get()
            result = prombsUtility(counts, self.data)
            self.queue_out.put((counts, result))
            self.queue_in.task_done()

def precomputeUtilityRec(counts, data, m, hashutil, queue, i_, j_):
    if m == 0:
        key = computeKey(counts)
        if not hashutil.get(key):
            queue.put(counts)
    elif m > 0:
        for i in range(i_, data['L']):
            for j in range(j_, data['K']):
                counts[j][i] += 1
                precomputeUtilityRec(counts, data, m-1, hashutil, queue, i, j)
                counts[j][i] -= 1

def precomputeUtility(counts, data, m, hashutil):
    queue_in  = Queue.Queue()
    queue_out = Queue.Queue()
    precomputeUtilityRec(counts, data, m, hashutil, queue_in, 0, 0)
    for i in range(options['threads']):
        t = ThreadUtility(queue_in, queue_out, data)
        t.setDaemon(True)
        t.start()
    queue_in.join()
    while not queue_out.empty():
        counts, result = queue_out.get()
        key = computeKey(counts)
        hashutil[key] = result

# precompute expectation
# ------------------------------------------------------------------------------

class ThreadExp(threading.Thread):
    def __init__(self, queue_in, queue_out, data):
        threading.Thread.__init__(self)
        self.queue_in  = queue_in
        self.queue_out = queue_out
        self.data      = data
    def run(self):
        while True:
            counts = self.queue_in.get()
            result = [ prombsExpectation(y, counts, self.data) for y in range(0, self.data['K']) ]
            self.queue_out.put((counts, result))
            self.queue_in.task_done()

def precomputeExpectationRec(counts, data, m, hashexp, queue, i_, j_):
    if m == 0:
        key = computeKey(counts)
        if not hashexp.get(key):
            queue.put(counts)
    elif m > 0:
        for i in range(i_, data['L']):
            for j in range(j_, data['K']):
                counts[j][i] += 1
                precomputeExpectationRec(counts, data, m-1, hashexp, queue, i, j)
                counts[j][i] -= 1

def precomputeExpectation(counts, data, m, hashexp):
    queue_in  = Queue.Queue()
    queue_out = Queue.Queue()
    precomputeExpectationRec(counts, data, m-1, hashexp, queue_in, 0, 0)
    for i in range(options['threads']):
        t = ThreadExp(queue_in, queue_out, data)
        t.setDaemon(True)
        t.start()
    queue_in.join()
    while not queue_out.empty():
        counts, result = queue_out.get()
        key = computeKey(counts)
        hashexp[key] = result

# main method to compute utilities
# ------------------------------------------------------------------------------

def computeUtility(counts, data):
    precomputeUtility(counts, data, options['look_ahead'], hashutil)
    precomputeExpectation(counts, data, options['look_ahead'], hashexp)
    return recursiveUtility(counts, data, options['look_ahead'], {}, hashutil, hashexp)

def selectItem(counts, data):
    # compute utility
    if options['strategy'] == 'uniform':
        utility = map(operator.neg, map(sum, zip(*counts)))
    elif options['strategy'] == 'uniform-random':
        utility = [ 0.0 for i in range(0, len(counts[0])) ]
    elif options['strategy'] == 'entropy':
        utility = computeUtility(counts, data)
    elif options['strategy'] == 'effective-counts':
        utility = computeUtility(counts, data)
    elif options['strategy'] == 'variance':
        utility = computeUtility(counts, data)
    else:
        raise IOError('Unknown strategy: '+options['strategy'])
    # filter utility
    if options['filter']:
        gainFilter = None
        exec options['filter']
        if not gainFilter is None:
            utility = gainFilter(utility, map(sum, zip(*counts)), result)
    return selectRandom(argmax(utility)), utility

# experiment
# ------------------------------------------------------------------------------

def experiment(index, data, result, msocket):
    if data['gt']:
        # load previous state of the random generator
        if result['states']:
            random.setstate(result['states'][index])
        # lapsing
        if options['lapsing'] >= random.uniform(0.0, 1.0):
            if 0.5 >= random.uniform(0.0, 1.0):
                return 0 # success
            else:
                return 1 # failure
        # real psychometric function
        sample = random.uniform(0.0, 1.0)
        # save state of the random generator
        if result['states']:
            result['states'][index] = random.getstate()
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
    if not result['counts']:
        result['counts'] = [ list(np.repeat(0, data['L'])),
                             list(np.repeat(0, data['L'])) ]
    for i in range(0, options['samples']):
        print >> sys.stderr, "Sampling... %.1f%%" % ((float(i)+1)/float(options['samples'])*100)
        index, utility = selectItem(result['counts'], data)
        event  = experiment(index, data, result, msocket)
        result['samples'].append(index)
        result['entropy'].append(prombsEntropy(result['counts'], data))
        result['counts'][event][index] += 1
    index, utility = selectItem(result['counts'], data)
    # update result
    bin_result = bin(result['counts'], data, options)
    bin_result['counts']  = result['counts']
    bin_result['samples'] = result['samples']
    bin_result['entropy'] = result['entropy']
    bin_result['states']  = result['states']
    bin_result['utility'] = utility
    if options['port']:
        close_msocket(msocket)
    return bin_result

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
        data['gt'] = config.readVector(config_parser, 'Ground Truth', 'gt', float)
        data['L'] = len(data['gt'])
        data['alpha'], data['beta'], data['gamma'] = \
            config.getParameters(config_parser, 'Ground Truth', os.path.dirname(config_file), data['K'], data['L'])
        config.readStrategy(config_parser, 'Ground Truth', options)
        result = loadResult()
        if not result['states']:
            result['states'] = config.readSeeds(config_parser, 'Ground Truth', 'seeds')
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
    'samples'              : 0,
    'look_ahead'           : 0,
    'epsilon'              : 0.00001,
    'n_moments'            : 2,
    'mgs_samples'          : (100,2000),
    'marginal'             : 0,
    'marginal_step'        : 0.01,
    'marginal_range'       : (0.0,1.0),
    'which'                : 0,
    'lapsing'              : 0.0,
    'threads'              : 1,
    'stacksize'            : 256*1024,
    'algorithm'            : 'prombs',
    'strategy'             : 'entropy',
    'port'                 : None,
    'filter'               : None,
    'visualization'        : None,
    'load'                 : None,
    'save'                 : None,
    'savefig'              : None,
    'verbose'              : False,
    'prombsTest'           : False,
    'compare'              : False,
    'bprob'                : False,
    'utility'              : False,
    'multibin_entropy'     : False,
    'differential_entropy' : False,
    'effective_counts'     : False,
    'model_posterior'      : True,
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=", "marginal", "marginal-range=",
                      "marginal-step=", "which=", "epsilon=", "moments", "look-ahead=",
                      "savefig=", "lapsing=", "port=", "threads=", "stacksize=",
                      "strategy=", "differential-entropy", "multibin-entropy",
                      "algorithm=", "samples=", "mgs-samples", "no-model-posterior"]
        opts, tail = getopt.getopt(sys.argv[1:], "mr:s:k:n:bhvt", longopts)
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
        if o == "--strategy":
            options["strategy"] = a
        if o == "--differential-entropy":
            options["differential_entropy"] = True
        if o == "--multibin-entropy":
            options["multibin_entropy"] = True
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
    if options["strategy"] == "effective-counts":
        options["effective_counts"] = True
    if (options["strategy"] == "entropy"    and 
        not options["differential_entropy"] and
        not options["multibin_entropy"]):
        options["differential_entropy"] = True
    if len(tail) != 1:
        usage()
        return 1
    parseConfig(tail[0])
    return 0

if __name__ == "__main__":
    sys.exit(main())
