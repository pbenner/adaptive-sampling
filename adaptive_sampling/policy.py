#! /usr/bin/env python

# Copyright (C) 2012 Philipp Benner
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

import Queue
import copy
import interface
import random
import statistics
import sys
import threading

# call the interface
################################################################################

def utility(counts_v, data, bin_options):
    """Call the binning library."""
    events        = len(counts_v)
    counts        = statistics.countStatistic(counts_v)
    alpha         = data['alpha']
    beta          = data['beta']
    gamma         = data['gamma']
    return interface.utility(events, counts, alpha, beta, gamma, bin_options)

def utilityAt(i, counts_v, data, bin_options):
    """Call the binning library."""
    events        = len(counts_v)
    counts        = statistics.countStatistic(counts_v)
    alpha         = data['alpha']
    beta          = data['beta']
    gamma         = data['gamma']
    return interface.utilityAt(i, events, counts, alpha, beta, gamma, bin_options)

# tools
################################################################################

def computeKey(position, counts):
    return tuple([position]+map(tuple, counts))

# determine the value of a sampling path
################################################################################

def value(path, counts, data, bin_options, hashutil):
    if len(path) == 0:
        return 0.0

    # if necessary compute the local utility for this count statistic
    key = computeKey(path[0], counts)
    if not hashutil.get(key):
        hashutil[key] = utilityAt(path[0], counts, data, bin_options)
    # get the local utility from the hashmap
    (expectation, utility) = hashutil.get(key)

    for i in range(data['K']):
        counts[i][path[0]] += 1
        utility += expectation[i]*value(path[1:], counts, data, bin_options, hashutil)
        counts[i][path[0]] -= 1

    return utility

# optimize a sampling path similar to the policy iteration algorithm
################################################################################

def optimize_entry(i, path_value, path, counts, data, bin_options, hashutil):
    changed = False
    stimuli = range(len(counts[0]))
    stimuli.remove(path[i])
    path_prime = copy.deepcopy(path)
    for x in stimuli:
        path_prime[i]    = x
        path_value_prime = value(path_prime, counts, data, bin_options, hashutil)
        if path_value_prime > path_value:
            changed    = True
            path_value = path_value_prime
            path[i]    = x
    return (path_value, path, changed)

def optimize(path, counts, data, bin_options, hashutil, full=False):
    changed    = True
    path_value = value(path, counts, data, bin_options, hashutil)
    decisions  = range(len(path))
    if not full:
        decisions.remove(0)
    while changed:
        for i in decisions:
            (path_value, path, changed) = optimize_entry(i, path_value, path, counts, data, bin_options, hashutil)
    return (path_value, path)

def u_star(length, counts, data, bin_options):
    if length <= 1:
        return utility(counts, data, bin_options)[1]
    stimuli  = range(len(counts[0]))
    path     = [ random.choice(stimuli) for i in range(length) ]
    utility  = [ 0.0 for i in stimuli ]
    hashutil = {}
    for x in stimuli:
        path[0] = x
        (path_value, path) = optimize(path, counts, data, bin_options, hashutil)
        utility[x] = path_value
    return utility

# threaded optimization of sampling paths
################################################################################

class OptimizationThread(threading.Thread):
    def __init__(self, length, counts, data, bin_options, queue_in, queue_out):
        threading.Thread.__init__(self)
        self.length      = length
        self.counts      = copy.deepcopy(counts)
        self.data        = copy.deepcopy(data)
        self.bin_options = bin_options
        self.queue_in    = queue_in
        self.queue_out   = queue_out
        self.hashutil    = {}
    def run(self):
        stimuli = range(len(self.counts[0]))
        path    = [ random.choice(stimuli) for i in range(self.length) ]
        while True:
            # get stimulus from queue
            x = self.queue_in.get()
            if self.bin_options['verbose']:
                sys.stderr.write('Processing stimulus ' + str(x) + '.\n')
            # set first element of the path to this stimulus
            path[0] = x
            # optimize all other elements of the path
            (path_value, path) = optimize(path, self.counts, self.data, self.bin_options, self.hashutil)
            # push result
            self.queue_out.put((x, path_value))
            self.queue_in.task_done()

def threaded_u_star(length, counts, data, bin_options):
    if length <= 1:
        return utility(counts, data, bin_options)[1]
    utility_queue_in  = Queue.Queue()
    utility_queue_out = Queue.Queue()
    utility_threads   = []
    stimuli = range(len(counts[0]))
    utility_vector = [ 0.0 for i in stimuli ]
    # launch daemon threads
    for i in range(bin_options['threads']):
        t = OptimizationThread(length, counts, data, bin_options,
                               utility_queue_in, utility_queue_out)
        t.setDaemon(True)
        t.start()
        utility_threads += [t]
    # fill queue and start computation
    stimuli = range(len(counts[0]))
    for x in stimuli:
        utility_queue_in.put(x)
    # wait for threads
    utility_queue_in.join()
    # process results
    while not utility_queue_out.empty():
        x, path_value = utility_queue_out.get()
        utility_vector[x] = path_value
        utility_queue_out.task_done()
    return utility_vector

# test functions
################################################################################

def test1(counts, data, bin_options):
    hashutil = {}

    for i in range(data['L']):
        print [i],
        print ": ",
        print value([i], counts, data, bin_options, hashutil)

    for i in range(data['L']):
        for j in range(data['L']):
            print [i, j],
            print ": ",
            print value([i, j], counts, data, bin_options, hashutil)

    for i in range(data['L']):
        for j in range(data['L']):
            for k in range(data['L']):
                print [i, j, k],
                print ": ",
                print value([i, j, k], counts, data, bin_options, hashutil)

    for i in range(data['L']):
        for j in range(data['L']):
            for k in range(data['L']):
                for l in range(data['L']):
                    print [i, j, k, l],
                    print ": ",
                    print value([i, j, k, l], counts, data, bin_options, hashutil)
