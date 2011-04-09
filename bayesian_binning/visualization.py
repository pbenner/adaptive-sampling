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

import math
import numpy as np
from matplotlib import *
from matplotlib.pyplot import *
from matplotlib.image import NonUniformImage
import matplotlib.patches as patches
import matplotlib.path as path

import statistics

font = {'family'     : 'serif',
        'weight'     : 'normal',
        'size'       : 12 }

smallfont = {'family'     : 'serif',
             'weight'     : 'normal',
             'size'       : 10 }

def plotGroundTruth(ax, x, gt):
    p = ax.plot(x, gt)
    ax.set_ylim(0,1)
    return p

def plotUtility(ax, x, result):
    p = ax.plot(x, result['differential_gain'])
    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    return p

def plotModelPosterior(ax, result):
    data = result['mpost'][:]
    N = len(data)
    x = np.arange(0, N+1, 1)
    data.insert(0, 0)
    ax.step(x, data, 'r--', where='mid', linewidth=1)
    ax.grid(True)
    left    = np.array(x[:-1]) + 0.5
    right   = np.array(x[1:])  + 0.5
    bottom  = np.zeros(len(left))
    top     = bottom + data[1:]
    XY      = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    barpath = path.Path.make_compound_path_from_polys(XY)
    patch   = patches.PathPatch(barpath, facecolor='green', edgecolor='gray', alpha=0.9)
    ax.add_patch(patch)
    ax.set_xlim(1, N)

def plotEffectiveCounts(ax, xorig, result):
    data = result['effective_counts'][:]
    x = list(xorig[:])
    x.insert(0, x[0]-1)
    data.insert(0, 0)
    N = len(data)
    p = ax.step(x, data, 'y--', where='mid', linewidth=1)
    ax.grid(True)
    left    = np.array(x[:-1]) + 0.5
    right   = np.array(x[1:])  + 0.5
    bottom  = np.zeros(len(left))
    top     = bottom + data[1:]
    XY      = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    barpath = path.Path.make_compound_path_from_polys(XY)
    patch   = patches.PathPatch(barpath, facecolor='yellow', edgecolor='gray', alpha=0.3)
    ax.add_patch(patch)
    ax.set_ylim(0, ax.get_ylim()[1]+1)
    return p

def plotCounts(ax, x, result):
    N = len(result['moments'][0])
    if len(result['samples']) > 0:
        n, bins, patches = ax.hist(result['samples'], N, range=(x[0]-0.5,x[-1]+0.5), normed=0, facecolor='yellow', alpha=0.8)
        ax.set_xlim(0, N-1)
        ax.set_ylim(0, ax.get_ylim()[1]+1)

def plotMultibinEntropy(ax, result):
    N = len(result['multibin_entropy'])
    x = np.arange(0, N+1, 1)
    result['multibin_entropy'].append(0)
    ax.step(x, result['multibin_entropy'], 'b--', where='mid', linewidth=1)

def plotSpikes(ax, x, timings):
    """Plot trials of spike trains."""
    X = []
    Y = []
    for i in range(0, len(timings)):
        for val in timings[i]:
            X.append(val)
            Y.append(i)
    ax.set_xlim(x[0],x[-1])
    ax.plot(X, Y, 'k|')

def plotBinBoundaries(ax, x, result):
    ax.plot(x[1:-1], result['bprob'][1:-1], 'g')
    ax.set_ylim(0,1)

def plotBin(ax, x, result):
    """Plot the binning result."""
    N = len(result['moments'][0])
    stddev = map(math.sqrt, statistics.centralMoments(result['moments'], 2))
    skew   = statistics.standardizedMoments(result['moments'], 3)
    ax.plot(x, [ a + b for a, b in zip(result['moments'][0], stddev) ], 'k--')
    ax.plot(x, [ a - b for a, b in zip(result['moments'][0], stddev) ], 'k--')
    p = ax.plot(x, result['moments'][0], 'r')
    ax.set_xlim(x[0],x[-1])
    [y1,y2] = ax.get_ylim()
    if y1 < 0: y1 = 0
    if y2 > 1: y2 = 1
    ax.set_ylim(y1, y2)
    return p

def plotMarginal(ax, x, result):
    N = len(result['moments'][0])
    y = np.linspace(0, 1, len(result['marginals'][0]))
    z = zip(*result['marginals'])
    im = NonUniformImage(ax, interpolation='nearest', cmap=cm.Greys)
    im.set_data(x, y, z)
    ax.images.append(im)

# ------------------------------------------------------------------------------

def plotBinning(result, options):
    x = np.arange(0, len(result['moments'][0]), 1)
    title = ''
    preplot  = None
    postplot = None
    if options['script']:
        exec options['script']
        if not preplot is None:
            x, title = preplot(result, options)
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    ax11 = fig.add_subplot(2,1,1)
    ax21 = fig.add_subplot(2,1,2)
    ax12 = ax11.twinx()
    ax22 = ax21.twinx()
    p11 = plotBin(ax11, x, result)
    p21 = plotModelPosterior(ax21, result)
    p12 = None
    p22 = None
    if result['marginals'] and options['marginal']:
        plotMarginal(ax11, x, result)
    if result['multibin_entropy'] and options['multibin_entropy']:
        plotMultibinEntropy(ax22, result)
    if result['bprob'] and options['bprob']:
        plotBinBoundaries(ax12, x, result)
    if result['differential_gain'] and options['differential_gain']:
        plotUtility(ax12, x, result)
    if result['effective_counts'] and options['effective_counts']:
        plotEffectiveCounts(ax22, x, result)
    if options['script'] and not postplot is None:
        postplot([ax11, ax12, ax21, ax22], [p11, p12, p21, p22],
                 result, options)

def plotBinningSpikes(x, timings, result, options):
    title = ''
    preplot  = None
    postplot = None
    if options['script']:
        exec options['script']
        if not preplot is None:
            title = preplot(result, options)
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    ax11 = fig.add_subplot(3,1,1, title=title)
    ax21 = fig.add_subplot(3,1,2)
    ax31 = fig.add_subplot(3,1,3)
    ax12 = ax11.twinx()
    ax22 = ax21.twinx()
    ax32 = ax31.twinx()
    p11 = plotSpikes(ax11, x, timings)
    p21 = plotBin   (ax21, x, result)
    p31 = plotModelPosterior(ax31, result)
    p12 = None
    p22 = None
    p32 = None
    if result['marginals'] and options['marginal']:
        plotMarginal(ax21, x, result)
    if result['multibin_entropy'] and options['multibin_entropy']:
        p32 = plotMutibinEntropy(ax32, result)
    if result['bprob'] and options['bprob']:
        p12 = plotBinBoundaries(ax12, x, result)
    if result['differential_gain'] and options['differential_gain']:
        p12 = plotUtility(ax12, result)
    if options['script'] and not postplot is None:
        postplot([ax11, ax12, ax21, ax22, ax31, ax32],
                 [p11, p12, p21, p22, p31, p31],
                 result, options)

def plotSampling(result, options, data):
    x = np.arange(0, data['bins'], 1)
    title = ''
    preplot  = None
    postplot = None
    if options['script']:
        exec options['script']
        if not preplot is None:
            x, title = preplot(result, options)
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    ax11 = fig.add_subplot(3,1,1, title=title)
    ax21 = fig.add_subplot(3,1,2)
    ax31 = fig.add_subplot(3,1,3)
    ax12 = ax11.twinx()
    ax22 = ax21.twinx()
    ax32 = ax31.twinx()
    p11 = plotBin(ax11, x, result)
    p21 = plotCounts(ax21, x, result)
    p31 = plotModelPosterior(ax31, result)
    p12 = None
    p22 = None
    p32 = None
    if data['gt']:
        p12 = plotGroundTruth(ax12, x, data['gt'])
    if result['marginals'] and options['marginal']:
        plotMarginal(ax11, x, result)
    if result['multibin_entropy'] and options['multibin_entropy']:
        plotMultibinEntropy(ax32, result)
    if result['bprob'] and options['bprob']:
        plotBinBoundaries(ax12, x, result)
    if result['differential_gain'] and options['strategy'] == 'differential-gain':
        p22 = plotUtility(ax22, x, result)
    if result['effective_counts'] and options['strategy'] == 'effective-counts':
        p22 = plotEffectiveCounts(ax22, x, result)
    if options['script'] and not postplot is None:
        postplot([ax11, ax12, ax21, ax22, ax31, ax32],
                 [p11, p12, p21, p22, p31, p31],
                 result, options)

def normalizeUtility(result):
    utility_matrix = []
    for utility in result['utility']:
        utility_min = min(utility)
        utility_max = max(utility)
        if utility_max - utility_min > 0 :
            utility_matrix.append(map(lambda x: (x-utility_min)/(utility_max-utility_min),utility))
        else:
            utility_matrix.append(map(lambda x: 0.5, utility))
    return utility_matrix

def plotUtilitySeries(result, options):
    utilitypreplot  = None
    utilitypostplot = None
    utility_matrix  = normalizeUtility(result)
    title = ''
    x = np.array(range(1, len(utility_matrix)+1))
    y = np.array(range(0, len(utility_matrix[0])))
    if options['script']:
        exec options['script']
        if not utilitypreplot is None:
            x, y, title = utilitypreplot(result)
    fig = figure()
    ax = fig.add_subplot(111, title=title)
    z = zip(*utility_matrix)
    p1 = NonUniformImage(ax, interpolation='nearest')
    p1.set_data(x, y, z)
    ax.images.append(p1)
    cb = colorbar(p1)
    z = [ result['samples'][i] for i in x ]
    ax.scatter(x, z, marker='o', c='w', s=10, zorder=10, edgecolors='none')
    ax.set_xlim(min(x)+0.5, max(x)+0.5)
    ax.set_ylim(min(y)-0.5, max(y)+0.5)
    if options['script'] and not utilitypostplot is None:
        utilitypostplot(ax, p1, cb)
