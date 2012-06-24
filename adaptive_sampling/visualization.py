# Copyright (C) 2010, 2011 Philipp Benner
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
import operator

import statistics

try:
    from matplotlib import *
    from matplotlib.pyplot import *
    from matplotlib.image import NonUniformImage
    import matplotlib.patches as patches
    import matplotlib.path as path
    import matplotlib.ticker as ticker
except ImportError:
    print "Error: Couldn't load matplotlib."
    exit(1)
except RuntimeError:
    print "Error: No x11 available."
    exit(1)

font = {'family'     : 'serif',
        'weight'     : 'normal',
        'size'       : 12 }

smallfont = {'family'     : 'serif',
             'weight'     : 'normal',
             'size'       : 10 }

def plotGroundTruth(ax, x, gt):
    p = ax.plot(x, gt)
    ax.set_ylim(0,1)
    return p[0]

def plotUtility(ax, x, result):
    p = ax.plot(x, result['utility'])
    ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    return p[0]

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
    data = map(operator.neg, result['utility'][:])
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
    return p[0]

def plotCounts(ax, x, result):
    N = len(result['moments'][0])
    if len(result['samples']) > 0:
        n, bins, patches = ax.hist(result['samples'], N, range=(x[0]-0.5,x[-1]+0.5), normed=0, facecolor='yellow', alpha=0.8)
        ax.set_xlim(0, N-1)
        ax.set_ylim(0, ax.get_ylim()[1]+1)

def plotSpikes(ax, x, timings):
    """Plot trials of spike trains."""
    X = []
    Y = []
    for i in range(0, len(timings)):
        for val in timings[i]:
            X.append(val)
            Y.append(i)
    ax.plot(X, Y, 'k|')
    ax.set_xlim(x[0],x[-1])

def plotBinBoundaries(ax, x, result):
    ax.plot(x[1:-1], result['bprob'][1:-1], 'g')
    ax.set_ylim(0,1)
    ax.set_xlim(x[0],x[-1])

def plotMoments(ax, x, result):
    """Plot the binning result."""
    N = len(result['moments'][0])
    stddev = map(math.sqrt, statistics.centralMoments(result['moments'], 2))
    stddev_upper = [ a + b for a, b in zip(result['moments'][0], stddev) ]
    stddev_lower = [ a - b for a, b in zip(result['moments'][0], stddev) ]
    ax.plot(x, stddev_upper, 'k-')
    ax.plot(x, stddev_lower, 'k-')
    ax.fill_between(x, stddev_lower, stddev_upper, linewidth=0, facecolor='red', alpha=0.3)
    p = ax.plot(x, result['moments'][0], 'r', linewidth=2)
    ax.set_xlim(x[0],x[-1])
    [y1,y2] = ax.get_ylim()
    if y1 < 0: y1 = 0
    if y2 > 1: y2 = 1
    ax.set_ylim(y1, y2)
    return p[0]

def plotSpikeRate(ax, x, result, dt):
    """Plot the binning result."""
    N = len(result['moments'][0])
    stddev = map(math.sqrt, statistics.centralMoments(result['moments'], 2))
    stddev_upper = [ (a + b)/dt for a, b in zip(result['moments'][0], stddev) ]
    stddev_lower = [ (a - b)/dt for a, b in zip(result['moments'][0], stddev) ]
    ax.plot(x, stddev_upper, 'k-')
    ax.plot(x, stddev_lower, 'k-')
    ax.fill_between(x, stddev_lower, stddev_upper, linewidth=0, facecolor='red', alpha=0.3)
    p = ax.plot(x, map(lambda x: x/dt, result['moments'][0]), 'r', linewidth=2)
    ax.set_xlim(x[0],x[-1])
    [y1,y2] = ax.get_ylim()
    if y1 < 0: y1 = 0
    ax.set_ylim(y1, y2)
    return p[0]

def plotMarginal(ax, x, result):
    N = len(result['moments'][0])
    y = np.linspace(0, 1, len(result['marginals'][0]))
    z = np.log(1+np.array(zip(*result['marginals'])))
    im = NonUniformImage(ax, norm=Normalize(0,5,clip=True), interpolation='nearest', cmap=cm.Greys)
    im.set_data(x, y, z)
    ax.images.append(im)

# ------------------------------------------------------------------------------

def plotBinning(result, options):
    x = np.arange(0, len(result['moments'][0]), 1)
    title = ''
    preplot  = None
    postplot = None
    if options['visualization']:
        exec options['visualization']
        if not preplot is None:
            x, title = preplot(result, options)
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    if (result['mpost'] and options['model_posterior']):
        ax11 = fig.add_subplot(2,1,1)
        ax21 = fig.add_subplot(2,1,2)
        ax12 = ax11.twinx()
        ax22 = ax21.twinx()
    else:
        ax11 = fig.add_subplot(1,1,1)
        ax21 = None
        ax12 = ax11.twinx()
        ax22 = None
    p11 = plotMoments(ax11, x, result)
    p12 = None
    p22 = None
    if result['marginals'] and options['marginal']:
        plotMarginal(ax11, x, result)
    if result['bprob'] and options['bprob']:
        p12 = plotBinBoundaries(ax12, x, result)
    if result['mpost'] and options['model_posterior']:
        p21 = plotModelPosterior(ax21, result)
    if options['visualization'] and not postplot is None:
        postplot([ax11, ax12, ax21, ax22], [p11, p12, p21, p22],
                 result, options)

def plotBinningSpikes(x, timings, result, options):
    dt    = 1.0
    title = ''
    preplot  = None
    postplot = None
    if options['visualization']:
        exec options['visualization']
        if not preplot is None:
            x, dt, timings, title = preplot(x, timings, result, options)
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    if result['mpost'] and options['model_posterior']:
        ax11 = fig.add_subplot(3,1,1, title=title)
        ax21 = fig.add_subplot(3,1,2)
        ax31 = fig.add_subplot(3,1,3)
        ax12 = ax11.twinx()
        ax22 = ax21.twinx()
        ax32 = ax31.twinx()
    else:
        ax11 = fig.add_subplot(2,1,1, title=title)
        ax21 = fig.add_subplot(2,1,2)
        ax31 = None
        ax12 = ax11.twinx()
        ax22 = ax21.twinx()
        ax32 = None
    p11 = plotSpikes(ax11, x, timings)
    p21 = plotSpikeRate(ax21, x, result, dt)
    p31 = None
    p12 = None
    p22 = None
    p32 = None
    if result['marginals'] and options['marginal']:
        plotMarginal(ax21, x, result)
    if result['bprob'] and options['bprob']:
        p12 = plotBinBoundaries(ax12, x, result)
    if result['mpost'] and options['model_posterior']:
        p31 = plotModelPosterior(ax31, result)
    if options['visualization'] and not postplot is None:
        postplot([ax11, ax12, ax21, ax22, ax31, ax32],
                 [p11, p12, p21, p22, p31, p31],
                 result, options)

def plotSampling(result, options, data):
    x = np.arange(0, data['L'], 1)
    title = ''
    preplot  = None
    postplot = None
    if options['visualization']:
        exec options['visualization']
        if not preplot is None:
            x, title = preplot(result, options)
    fig = figure(1, (8.0, 6.0))
    fig.subplots_adjust(left=0.09, bottom=0.06, right=0.89, top=0.96, hspace=0.35)
    if result['mpost'] and options['model_posterior']:
        ax11 = fig.add_subplot(3,1,1, title=title)
        ax21 = fig.add_subplot(3,1,2)
        ax31 = fig.add_subplot(3,1,3)
        ax12 = ax11.twinx()
        ax22 = ax21.twinx()
        ax32 = ax31.twinx()
    else:
        ax11 = fig.add_subplot(2,1,1, title=title)
        ax21 = fig.add_subplot(2,1,2)
        ax31 = None
        ax12 = ax11.twinx()
        ax22 = ax21.twinx()
        ax32 = None
    p11 = plotMoments(ax11, x, result)
    p21 = plotCounts(ax21, x, result)
    p31 = None
    p12 = None
    p22 = None
    p32 = None
    if result['marginals'] and options['marginal']:
        plotMarginal(ax11, x, result)
    if data['gt']:
        p12 = plotGroundTruth(ax12, x, data['gt'])
    if result['bprob'] and options['bprob']:
        p12 = plotBinBoundaries(ax12, x, result)
    if result['utility']:
        p22 = plotUtility(ax22, x, result)
        # if options['strategy'] == 'effective-counts':
        #     p22 = plotEffectiveCounts(ax22, x, result)
        # else:
        #     p22 = plotUtility(ax22, x, result)
    if result['mpost'] and options['model_posterior']:
        p31 = plotModelPosterior(ax31, result)
    if options['visualization'] and not postplot is None:
        postplot([ax11, ax12, ax21, ax22, ax31, ax32],
                 [p11, p12, p21, p22, p31, p31],
                 result, options)
