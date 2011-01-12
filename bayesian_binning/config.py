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

import numpy as np
import ConfigParser

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

def readAlpha(config_parser, n, section, converter):
    alpha = []
    if config_parser.has_option(section, 'alpha'):
        alpha = readVector(config_parser, section, 'alpha', converter)
        if len(alpha) < n:
            raise ValueError("Not enough alpha parameters.")
        if len(alpha) > n:
            raise ValueError("Too many alpha parameters.")
    else:
        alpha = list(np.repeat(1, n))
    return alpha

def readModelPrior(config_parser, n, section, converter):
    models = []
    mprior = list(np.repeat(1, n))
    if config_parser.has_option(section, 'mprior'):
        mprior = list(np.repeat(0, n))
        models = readVector(config_parser, section, 'mprior', converter)
        num_models = len(models)
        for model in models:
            mprior[model-1] = 1.0/num_models
    return mprior
