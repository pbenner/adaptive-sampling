#! /usr/bin/env python

# Copyright (C) 2010, 2011, 2012 Philipp Benner
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
import copy

from itertools import izip

def importMatplotlib(backend=None):
    global vis
    global use
    global pyplot
    if not globals().has_key('use'):
        from matplotlib import use
        if backend:
            use(backend)
    from matplotlib import pyplot
    if not globals().has_key('vis'):
        import visualization as vis

import config
import interface
import statistics
import policy

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
    print "       --distances                    - compute Kullback-Leibler distance for each sample"
    print "       --lapsing=p                    - specify a lapsing probability"
    print "       --look-ahead=N                 - recursion depth for the sampling look ahead"
    print "       --hmm                          - use hidden Markov model"
    print "       --rho                          - rho parameter for the HMM"
    print "   -m  --density                     - compute full density distribution"
    print "   -r  --density-range=(FROM,TO)     - limit range for the density distribution"
    print "   -s  --density-step=STEP           - step size for the density distribution"
    print "       --no-model-posterior           - do not compute the model posterior"
    print "       --epsilon=EPSILON              - epsilon for the extended prombs"
    print "   -n  --samples=N                    - number of samples"
    print "   -k  --moments=N                    - compute the first N>=2 moments"
    print "       --which=EVENT                  - for which event to compute the binning"
    print "       --algorithm=NAME               - select an algorithm [mgs, default: prombs]"
    print "       --mgs-samples=BURN_IN:SAMPLES  - number of samples [default: 100:2000] for mgs"
    print "       --path-iteratin                - use path iteration algorithm instead of backward"
    print "                                        to compute n-step utilities"
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
    print "       --video=BASENAME               - generate a video from the sampling process"
    print
    print "   -h, --help                         - print help"
    print "   -v, --verbose                      - be verbose"
    print "   -t, --prombsTest                   - test prombs algorithm"
    print
    print "Sampling strategies:"
    print "       --strategy=STRATEGY            - uniform, uniform-random, kl-divergence (default),"
    print "                                        effective-counts, effective-posterior-counts"
    print "       --kl-psi                       - if kl-divergence is selected as strategy then use"
    print "                                        add the psi divergence"
    print "       --kl-multibin                  - if kl-divergence is selected as strategy then use"
    print "                                        add the multibin divergence"
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

def roundArray(array, precision):
    return map(lambda x: round(x, precision), array)

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
    config.set('Sampling Result', 'density',   "\n"+"\n".join(map(lambda arg: " ".join(map(str, arg)), result['density'])))
    config.set('Sampling Result', 'samples',   " ".join(map(str, result['samples'])))
    config.set('Sampling Result', 'utility',   " ".join(map(str, result['utility'])))
    config.set('Sampling Result', 'bprob',     " ".join(map(str, result['bprob'])))
    config.set('Sampling Result', 'mpost',     " ".join(map(str, result['mpost'])))
    config.set('Sampling Result', 'distances', " ".join(map(str, result['distances'])))
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

        distances = []
        bprob     = []

        counts    = config.readMatrix(config_parser, 'Sampling Result', 'counts',  int)
        moments   = config.readMatrix(config_parser, 'Sampling Result', 'moments', float)
        density   = config.readMatrix(config_parser, 'Sampling Result', 'density', float)
        samples   = config.readVector(config_parser, 'Sampling Result', 'samples', int)
        mpost     = config.readVector(config_parser, 'Sampling Result', 'mpost',   float)
        states    = config.readStates(config_parser, 'Sampling Result', 'states')
        if config_parser.has_option('Sampling Result', 'distances'):
            distances = config.readVector(config_parser, 'Sampling Result', 'distances', float)
        if config_parser.has_option('Sampling Result', 'bprob'):
            bprob = config.readVector(config_parser, 'Sampling Result', 'bprob',   float)
        result = {
            'distances' : distances,
            'moments'   : moments,
            'density'   : density,
            'bprob'     : bprob,
            'mpost'     : mpost,
            'counts'    : counts,
            'samples'   : samples,
            'states'    : states }
    else:
        result = {
            'distances' : [],
            'moments'   : [],
            'density'   : [],
            'bprob'     : [],
            'mpost'     : [],
            'counts'    : [],
            'samples'   : [],
            'states'    : [] }
    return result

# call the interface
# ------------------------------------------------------------------------------

def call_posterior(counts_v, data, bin_options):
    """Call the binning library."""
    events = len(counts_v)
    counts = statistics.countStatistic(counts_v)
    alpha  = data['alpha']
    beta   = data['beta']
    gamma  = data['gamma']
    return interface.posterior(events, counts, alpha, beta, gamma, bin_options)

def call_utility(counts_v, data, bin_options):
    """Call the binning library."""
    events        = len(counts_v)
    counts        = statistics.countStatistic(counts_v)
    alpha         = data['alpha']
    beta          = data['beta']
    gamma         = data['gamma']
    return interface.utility(events, counts, alpha, beta, gamma, bin_options)

def call_distance(x, y, counts_v, data, bin_options):
    """Call the binning library."""
    events        = len(counts_v)
    counts        = statistics.countStatistic(counts_v)
    alpha         = data['alpha']
    beta          = data['beta']
    gamma         = data['gamma']
    return interface.distance(x, y, events, counts, alpha, beta, gamma, bin_options)

# recursive computation of utilities
# ------------------------------------------------------------------------------

def computeKey(counts):
    return tuple(map(tuple, counts))

def recursiveUtility(counts, data, m, hashutil):
    key = computeKey(counts)
    # compute immediate utility if not done yet
    if not hashutil.get(key):
        hashutil[key] = call_utility(counts, data, options)
    # receive utility and expectation
    tmp = hashutil[key]
    (expectation, utility) = (tmp['expectation'], tmp['utility'])
    if m == 0:
        return utility
    else:
        result = [ 0.0 ] * data['L']
        # for all stimuli
        for x in range(0, data['L']):
            # loop over events to compute the expectation
            for y in range(0, data['K']):
                counts[y][x] += 1
                result[x] += expectation[y][x]*max(recursiveUtility(counts, data, m-1, hashutil))
                counts[y][x] -= 1
        # add the utility for the result
        result = map(sum, zip(utility, result))
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
            result = call_utility(counts, self.data, options)
            self.queue_out.put((counts, result))
            self.queue_in.task_done()

def precomputeUtilityRec(counts, data, m, hashutil, queue, i_, j_):
    # put counts into the hashmap
    if m >= 0:
        key = computeKey(counts)
        if not hashutil.get(key):
            queue.put(copy.deepcopy(counts))
    # if recursion limit has not been reached then nest further
    if m > 0:
        for j in range(j_, data['K']):
            if j == j_:
                i_from = i_
            else:
                i_from = 0
            for i in range(i_from, data['L']):
                counts[j][i] += 1
                precomputeUtilityRec(counts, data, m-1, hashutil, queue, i, j)
                counts[j][i] -= 1

utility_threads   = []
utility_queue_in  = Queue.Queue()
utility_queue_out = Queue.Queue()
def precomputeUtility(counts, data, m, hashutil):
    global utility_threads
    if not utility_threads:
        # launch daemon threads
        for i in range(options['threads']):
            t = ThreadUtility(utility_queue_in, utility_queue_out, data)
            t.setDaemon(True)
            t.start()
            utility_threads += [t]
    precomputeUtilityRec(counts, data, m, hashutil, utility_queue_in, 0, 0)
    utility_queue_in.join()
    while not utility_queue_out.empty():
        counts, result = utility_queue_out.get()
        key = computeKey(counts)
        hashutil[key] = result
        utility_queue_out.task_done()

# main method to compute utilities
# ------------------------------------------------------------------------------

def computeUtility(counts, data):
    bin_options = options.copy()
    hashutil    = {}
    if options['path_iteration']:
        # print policy.optimize([1,1,1,1], counts, data, bin_options)
        # policy.test1(counts, data, bin_options)
        utility = policy.threaded_u_star(options['look_ahead']+1, counts, data, bin_options)
    else:
        precomputeUtility(counts, data, options['look_ahead'], hashutil)
        utility = recursiveUtility(counts, data, options['look_ahead'], hashutil)
    return utility

def selectItem(counts, data):
    # compute utility
    if options['strategy'] == 'uniform':
        utility = map(operator.neg, map(sum, zip(*counts)))
    elif options['strategy'] == 'uniform-random':
        utility = [ 0.0 for i in range(0, len(counts[0])) ]
    elif options['strategy'] == 'kl-divergence':
        utility = computeUtility(counts, data)
    elif options['strategy'] == 'kl-multibin':
        utility = computeUtility(counts, data)
    elif options['strategy'] == 'effective-counts':
        utility = computeUtility(counts, data)
    elif options['strategy'] == 'effective-posterior-counts':
        utility = computeUtility(counts, data)
    else:
        raise IOError('Unknown strategy: '+options['strategy'])
    # filter utility
    if options['filter']:
        gainFilter = None
        exec options['filter']
        if not gainFilter is None:
            utility = gainFilter(utility, map(sum, zip(*counts)))
#    utility = roundArray(utility, 10)
    if options['verbose']:
        sys.stderr.write(str(utility) + '.\n')
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

# video
# ------------------------------------------------------------------------------

def save_frame(result, data, utility, i):
    importMatplotlib('Agg')
    from matplotlib.pyplot import savefig
    pyplot.clf()
    bin_result = call_posterior(result['counts'], data, options)
    bin_result['counts']  = result['counts']
    bin_result['samples'] = result['samples']
    bin_result['states']  = result['states']
    bin_result['utility'] = utility
    vis.plotSampling(bin_result, options, data)
    savefig('%s_%03d.png' % (options['video'], i), bbox_inches='tight', pad_inches=0)

def save_video():
    import subprocess
    command = ('mencoder',
               'mf://%s_*.png' % options['video'],
               '-mf',
               'type=png:fps=4',
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               '%s.avi' % options['video'])
    subprocess.check_call(command)

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
        event = experiment(index, data, result, msocket)
        # compute Kullback-Leibler distance for the new sample
        if options['distances']:
            result['distances'].append(call_distance(index, event, result['counts'], data, options))
        # record new sample
        result['samples'  ].append(index)
        result['counts'   ][event][index] += 1
        if options['video']:
            save_frame(result, data, utility, i)
    index, utility = selectItem(result['counts'], data)
    # update result
    bin_result = call_posterior(result['counts'], data, options)
    bin_result['counts']    = result['counts']
    bin_result['distances'] = result['distances']
    bin_result['samples']   = result['samples']
    bin_result['states']    = result['states']
    bin_result['utility']   = utility
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
    if options['video']:
        save_video()
    if options['save']:
        saveResult(result)
    else:
        if options['savefig']:
            importMatplotlib('Agg')
            from matplotlib.pyplot import savefig
            vis.plotSampling(result, options, data)
#            savefig(options['savefig'], bbox_inches='tight', pad_inches=0)
            savefig(options['savefig'])
        else:
            importMatplotlib()
            from matplotlib.pyplot import show
            vis.plotSampling(result, options, data)
            show()

# main
# ------------------------------------------------------------------------------

options = {
    'samples'                    : 0,
    'look_ahead'                 : 0,
    'epsilon'                    : 0.00001,
    'n_moments'                  : 2,
    'mgs_samples'                : (100,2000),
    'density'                    : 0,
    'density_step'               : 0.01,
    'density_range'              : (0.0,1.0),
    'which'                      : 0,
    'lapsing'                    : 0.0,
    'threads'                    : 1,
    'stacksize'                  : 256*1024,
    'algorithm'                  : 'prombs',
    'strategy'                   : 'kl-divergence',
    'video'                      : None,
    'port'                       : None,
    'filter'                     : None,
    'visualization'              : None,
    'load'                       : None,
    'save'                       : None,
    'savefig'                    : None,
    'verbose'                    : False,
    'prombsTest'                 : False,
    'compare'                    : False,
    'bprob'                      : False,
    'utility'                    : False,
    'kl_psi'                     : False,
    'kl_multibin'                : False,
    'effective_counts'           : False,
    'effective_posterior_counts' : False,
    'model_posterior'            : True,
    'distances'                  : False,
    'hmm'                        : False,
    'path_iteration'             : False,
    'rho'                        : 0.4
    }

def main():
    global options
    try:
        longopts   = ["help", "verbose", "load=", "save=", "density", "density-range=",
                      "density-step=", "which=", "epsilon=", "moments", "look-ahead=",
                      "savefig=", "lapsing=", "port=", "threads=", "stacksize=",
                      "strategy=", "kl-psi", "kl-multibin", "algorithm=", "samples=",
                      "mgs-samples", "no-model-posterior", "video=", "hmm", "rho=",
                      "path-iteration", "distances" ]
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
        if o in ("-m", "--density"):
            options["density"] = True
        if o in ("-r", "--density-range"):
            options['density_range'] = tuple(map(float, a[1:-1].split(',')))
        if o in ("-s", "--density-step"):
            options["density_step"] = float(a)
        if o in ("-k", "--moments"):
            if int(a) >= 2:
                options["n_moments"] = int(a)
            else:
                usage()
                return 0
        if o in ("-h", "--help"):
            usage()
            return 0
        if o == "--distances":
            options["distances"] = True
        if o == "-b":
            options["bprob"] = True
        if o == "--lapsing":
            options["lapsing"] = float(a)
        if o == "--strategy":
            options["strategy"] = a
        if o == "--kl-psi":
            options["kl_psi"] = True
        if o == "--kl-multibin":
            options["kl_multibin"] = True
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
        if o == "--video":
            options["video"] = a
        if o == "--hmm":
            options["hmm"] = True
        if o == "--rho":
            options["rho"] = float(a)
        if o == "--path-iteration":
            options["path_iteration"] = True
    if (options["strategy"] == "kl-divergence" and
        options["kl_psi"]   == False           and
        options["kl_multibin"]  == False):
       # set default kl divergence
       options["kl_psi"]   = True
    if options["strategy"] == "effective-counts":
        options["effective_counts"] = True
    if options["strategy"] == "effective-posterior-counts":
        options["effective_posterior_counts"] = True
    if len(tail) != 1:
        usage()
        return 1
    interface.init(options['epsilon'])
    parseConfig(tail[0])
    interface.free()
    return 0

if __name__ == "__main__":
    sys.exit(main())
