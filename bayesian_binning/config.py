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

import bayesian_binning.statistics    as statistics

## functions for reading config options
################################################################################

def readVector(config, section, option, converter):
    vector_str = config.get(section, option)
    vector     = map(converter, vector_str.split(' '))
    return vector

def readMatrix(config, section, option, converter):
    matrix_str = config.get(section, option)
    matrix     = []
    for line in matrix_str.split('\n'):
        if line != '':
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

## helper functions for setting prior parameters
################################################################################

def generate_alpha(alpha_v):
    K = len(alpha_v)
    L = len(alpha_v[0])
    alpha = np.zeros([K, L, L])
    for ks in range(0, L):
        for ke in range(ks, L):
            c = np.zeros(K)
            for i in range(ks, ke+1):
                for j in range(0, K):
                    c[j] += alpha_v[j][i] / float(ke - ks + 1)
            for j in range(0, K):
                alpha[j][ks][ke] = c[j]
    return alpha

def generate_beta(beta_v, num_models):
    if beta_v == []:
        return np.ones(num_models)/num_models
    beta = np.zeros(num_models)
    free_models = 0.0
    p_sum = 0.0
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
    return beta

def generate_gamma(num_models):
    gamma = np.ones([num_models, num_models])
    return np.triu(gamma)

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

def readFilter(config_parser, section, dir):
    return readScript(config_parser, section, dir, 'filter')

def readVisualization(config_parser, section, dir):
    return readScript(config_parser, section, dir, 'visualization')
