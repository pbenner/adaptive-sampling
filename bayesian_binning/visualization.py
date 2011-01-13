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
import matplotlib.patches as patches
import matplotlib.path as path

import statistics

font = {'family'     : 'serif',
        'color'      : 'k',
        'weight'     : 'normal',
        'size'       : 12 }

def plotGroundTruth(ax, gt):
    x = np.arange(0, len(gt), 1)
    ax.plot(x, gt, label='Ground Truth')
    ax.set_ylabel('Ground Truth', font)
    ax.legend(loc=2, frameon=False)

def plotUtility(ax, result):
    x = np.arange(0, len(result['differential_gain']), 1)
    plot(x, result['differential_gain'], label=r'$U_x$')
    ax.set_ylabel(r'$U_x$', font)
    ax.ticklabel_format(style='sci', scilimits=(0,0))
    ax.legend(loc=4, frameon=False)

def plotModelPosterior(ax, result):
    N = len(result['mpost'])
    x = np.arange(0, N+1, 1)
    result['mpost'].insert(0, 0)
    ax.step(x, result['mpost'], 'r--', where='mid', linewidth=1)
    ax.grid(True)

    left    = np.array(x[:-1]) + 0.5
    right   = np.array(x[1:])  + 0.5
    bottom  = np.zeros(len(left))
    top     = bottom + result['mpost'][1:]
    XY      = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    barpath = path.Path.make_compound_path_from_polys(XY)
    patch   = patches.PathPatch(barpath, facecolor='green', edgecolor='gray', alpha=0.8)
    ax.add_patch(patch)
    ax.set_xlabel(r'$m_B$',  font)
    ax.set_ylabel(r'$P(m_B|E)$', font)
    ax.set_xlim(0, N)

def plotCounts(ax, result):
    N = len(result['mpost'])
    if len(result['samples']) > 0:
        n, bins, patches = ax.hist(result['samples'], N, normed=0, facecolor='yellow', alpha=0.8)
        ax.set_ylim(0, ax.get_ylim()[1]+1)
    ax.set_ylabel('Counts', font)
    
def plotMultibinEntropy(ax, result):
    N = len(result['multibin_entropy'])
    x = np.arange(0, N+1, 1)
    result['multibin_entropy'].append(0)
    ax.step(x, result['multibin_entropy'], 'b--', where='mid', linewidth=1)
    ax.set_ylabel(r'$H(\mathcal{B}|E,m_B)$', font)

def plotSpikes(ax, x, timings):
    """Plot trials of spike trains."""
    X = []
    Y = []
    for i in range(0, len(timings)):
        for val in timings[i]:
            X.append(val)
            Y.append(i)
    ax.set_xlim(x[0],x[-1])
    ax.set_ylabel(r'$t$')
    ax.set_ylabel('Trial')
    ax.plot(X, Y, 'k|')

def plotBinBoundaries(ax, x, result):
    ax.plot(x[1:-1], result['bprob'][1:-1], 'g')
    ax.set_ylim(0,1)
    ax.set_ylabel(r'$P(\Rsh_i|E)$', font)
#    nbins = argmax(result['mpost'])[1]+1
#    bprob_max = sorted(result['bprob'])[-nbins:]
#    for b in bprob_max:
#        i = bprob.index(b)
#        ax.axvline(x[i])

def plotBin(ax, x, result):
    """Plot the binning result."""
    N = len(result['moments'][0])
    if x==None:
        x = np.arange(0, N, 1)
    ax.set_xlim(x[0],x[-1])
    ax.set_ylim(0,1)
    stddev = map(math.sqrt, statistics.centralMoments(result['moments'], 2))
    skew     = statistics.standardizedMoments(result['moments'], 3)
#    kurtosis = statistics.standardizedMoments(result['moments'], 4)
#    tail     = statistics.standardizedMoments(result['moments'], 5)
    ax.plot(x, [ a + b for a, b in zip(result['moments'][0], stddev) ], 'k--')
    ax.plot(x, [ a - b for a, b in zip(result['moments'][0], stddev) ], 'k--')
#    ax.plot(x, [ a - b for a, b in zip(result['moments'][0], skew) ], 'g--')
#    ax.plot(x, [ a + b for a, b in zip(result['moments'][0], tail) ], 'y--')
    ax.plot(x, result['moments'][0], 'r', label='$P(S_i|E)$')
    ax.set_xlabel('t',  font)
    ax.set_ylabel(r'$P(S_i|E)$', font)
    ax.legend(loc=4, frameon=False)

def plotBinning(x, result, plot_bprob, plot_multientropy):
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    plotBin(ax1, None, result)
    plotModelPosterior(ax2, result)
    if result['multibin_entropy'] and plot_multientropy:
        plotMultibinEntropy(ax2.twinx(), result)
    if result['bprob'] and plot_bprob:
        plotBinBoundaries(ax1.twinx(), x, result)
    if result['differential_gain']:
        plotUtility(ax1.twinx(), result)
    show()

def plotSampling(x, result, gt, options):
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    ax1 = fig.add_subplot(3,1,1)
    ax2 = fig.add_subplot(3,1,2)
    ax3 = fig.add_subplot(3,1,3)
    title(str(len(result['samples']))+" Samples")
    plotBin(ax1, None, result)
    plotGroundTruth(ax1.twinx(), gt)
    plotCounts(ax2, result)
    plotModelPosterior(ax3, result)
    if result['multibin_entropy'] and options['multibin_entropy']:
        plotMultibinEntropy(ax3.twinx(), result)
    if result['bprob'] and options['bprob']:
        plotBinBoundaries(ax1.twinx(), x, result)
    if result['differential_gain'] and options['differential_gain']:
        plotUtility(ax2.twinx(), result)
    show()

def plotBinningSpikes(x, timings, result, options):
    fig = figure()
    fig.subplots_adjust(hspace=0.35)
    ax1 = fig.add_subplot(3,1,1)
    ax2 = fig.add_subplot(3,1,2)
    ax3 = fig.add_subplot(3,1,3)
    plotSpikes(ax1, x, timings)
    plotBin   (ax2, x, result)
    plotModelPosterior(ax3, result)
    if result['multibin_entropy'] and options['multibin_entropy']:
        plotMutibinEntropy(ax3.twinx(), result)
    if result['bprob'] and options['bprob']:
        plotBinBoundaries(ax1.twinx(), x, result)
    if result['differential_gain'] and options['differential_gain']:
        plotUtility(ax1.twinx(), result)
    show()
