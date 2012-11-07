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

import os.path
import numpy as np
import ConfigParser
import re
import random
import math

import statistics

from interface import get_huge_val

## functions for reading config options
################################################################################

def readVector(config, section, option, converter):
    vector_str = config.get(section, option)
    vector_str.strip()
    if vector_str != '':
        vector = map(converter, vector_str.split(' '))
    else:
        vector = []
    return vector

def readMatrix(config, section, option, converter):
    matrix_str = config.get(section, option)
    matrix     = []
    for line in matrix_str.split('\n'):
        if line != '':
            line.strip()
            matrix.append([converter(a) for a in line.split(' ')])
    return matrix

def readStrategy(config_parser, section, options):
    if config_parser.has_option(section, 'strategy'):
        options['strategy'] = config_parser.get(section, 'strategy')

def readAlpha(config_parser, events, bins, section, converter):
    alpha = []
    if config_parser.has_option(section, 'alpha'):
        alpha = readMatrix(config_parser, section, 'alpha', converter)
        if not len(alpha) == events or not len(alpha[0]) == bins:
            raise ValueError("Wrong number of alpha parameters.")
    else:
        alpha = np.repeat(1, events*bins).reshape(events,bins)
    return alpha

def readCounts(config_parser, section):
    counts = readMatrix(config_parser, section, 'counts', float)
    return statistics.countStatistic(counts)

def readSeeds(config_parser, section, option):
    if config_parser.has_option(section, option):
        seeds = readVector(config_parser, section, option, float)
        rand_states = []
        for i in seeds:
            random.seed(i)
            rand_states.append(random.getstate())
    else:
        rand_states = []
    return rand_states

def readStates(config_parser, section, option):
    states = []
    if config_parser.has_option(section, option):
        states_str = config_parser.get(section, option)
        for state_str in states_str.split('\n'):
            if state_str != '':
                states.append(eval(state_str))
    return states

## helper functions for setting prior parameters
################################################################################

def generate_alpha(alpha_v):
    """Generate a set of default alpha parameters."""
    K = len(alpha_v)
    L = len(alpha_v[0])
    ones  = np.ones(L, dtype=float)
    alpha = np.zeros([K, L, L])
    for k in range(0, K):
        alpha[k] = np.triu(np.outer(ones,alpha_v[k])).cumsum(axis=1)
        for i in range(L):
            for j in range(i,L):
                alpha[k][i,j] /= float(j-i+1)
    return alpha

def generate_beta(beta_v, num_models, apply_model_penalty=True):
    """Generate a set of default beta parameters on log scale."""
    if beta_v == []:
        beta = np.ones(num_models)/num_models
    else:
        beta = np.zeros(num_models)
    free_models = 0.0
    p_sum       = 0.0
    for elem in beta_v:
        if isinstance(elem, tuple):
            p_sum += elem[1]
        else:
            free_models += 1.0
    for elem in beta_v:
        if isinstance(elem, tuple):
            beta[elem[0]-1] = float(elem[1])
        else:
            beta[elem-1] = (1.0-p_sum)/free_models
    # convert to log scale
    for i in range(num_models):
        if beta[i] > 0.0:
            beta[i] = math.log(beta[i])
        else:
            beta[i] = -get_huge_val()
    # multiply beta with (L-1 over m-1)^-1 on log scale
    if apply_model_penalty:
        model_penalty = [ statistics.gsl_sf_lnchoose(num_models-1, m) for m in range(num_models) ]
        beta = map(lambda (a, b): a - b, zip(beta, model_penalty))
    return beta

def generate_gamma(num_models):
    """Generate a set of default gamma parameters."""
    gamma = np.ones([num_models, num_models])
    return np.triu(gamma)

def generate_counts(events):
    K = len(events)
    L = len(events[0])
    ones   = np.ones(L, dtype=float)
    counts = np.zeros([K, L, L])
    for k in range(0, K):
        counts[k] = np.triu(np.outer(ones,events[k])).cumsum(axis=1)
    return counts

## read additional scripts
################################################################################

def readScript(config_parser, section, dir, name):
    if config_parser.has_option(section, name):
        file = config_parser.get(section, name)
        f = open(os.path.abspath(dir)+'/'+file, 'r')
        str = f.read()
        f.close()
        return str
    else:
        return None

def getParameters(config_parser, section, dir, K, L):
    script = readScript(config_parser, section, dir, 'parameters')
    if script:
        parameters = None
        exec script
        if not parameters is None:
            return parameters(K, L)
    else:
        alpha = generate_alpha(np.ones([K, L]))
        beta  = generate_beta([], L)
        gamma = generate_gamma(L)
        return alpha, beta, gamma

def readFilter(config_parser, section, dir, options):
    options['filter'] =  readScript(config_parser, section, dir, 'filter')

def readVisualization(config_parser, section, dir, options):
    options['visualization'] = readScript(config_parser, section, dir, 'visualization')

def readAlgorithm(config_parser, section, dir, options):
    if config_parser.has_option(section, 'algorithm'):
        options['algorithm'] = config_parser.get(section, 'algorithm')

def readMgsSamples(config_parser, section, dir, options):
    if config_parser.has_option(section, 'mgs-samples'):
         samples_str = config_parser.get(section, 'mgs-samples')
         options['mgs_samples'] = tuple(map(int, samples_str.split(":")))
