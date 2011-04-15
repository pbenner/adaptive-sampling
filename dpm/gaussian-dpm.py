#! /usr/bin/env python

# imports
################################################################################

import math
import time
import gobject
import gtk

import numpy               as     np
import numpy.random.mtrand as     mt
import numpy.random        as     rd

from   matplotlib          import *
use('GTKAgg')

from   matplotlib.pyplot   import *
from   matplotlib.image    import NonUniformImage
import matplotlib.patches  as     patches
import matplotlib.path     as     path

# local imports
################################################################################

from interface  import *

# Gaussian DPM
################################################################################

def biNormalDensity(mu, cov):
    return (lambda X, Y: mlab.bivariate_normal(X, Y, np.sqrt(cov[0,0]), np.sqrt(cov[1,1]), mu[0], mu[1], cov[0,1]))

class GaussianDPM():
    def __init__(self, n, k):
        dpm_init(n, k)
        self.cluster_colors = [ tuple(rd.rand(3)) for i in range(0, n*k) ]
        self.steps = 0
    def print_clusters(self):
        dpm_print()
    def num_clusters(self):
        return dpm_num_clusters()
    def sample(self, n):
        dpm_sample(n)
        self.steps += 1
    def hist_means():
        return dpm_hist_means()
    def means():
        return dpm_means()
    def plotData(self, ax):
        num_clusters = dpm_num_clusters()
        for c in range(0, num_clusters):
            x, y = zip(*dpm_cluster(c))
            original_tags = dpm_original_tags(c)
            for (x_, y_, c_) in zip(x, y, original_tags):
                ax.scatter(x_, y_, c=self.cluster_colors[c_])
    def plotResult(self, ax):
        num_clusters = dpm_num_clusters()
        for c in range(0, num_clusters):
            x, y = zip(*dpm_cluster(c))
            ax.scatter(x, y, c=self.cluster_colors[c])
            ax.set_title("Clustering Result K="+str(num_clusters)+", N="+str(self.steps))
    def plotJoint(self, ax):
        im = NonUniformImage(ax, interpolation='bilinear', cmap=cm.gray)
        im.set_data(self.x, self.y, self.Z)
        ax.images.append(im)
    def plotStatistics(self, ax1):
        ax2 = ax1.twinx()
        p1  = ax1.plot(dpm_hist_switches())
        p2  = ax2.plot(dpm_hist_likelihood(), color='green')
        ax1.set_ylabel("Class switches")
        ax2.set_ylabel("Likelihood")
        ax1.set_xlabel("iteration")
        ax1.legend([p1, p2], ["Mean class switches", "Mean likelihood"])

class InteractiveGDPM(GaussianDPM):
    def __init__(self, n, k, ax):
        GaussianDPM.__init__(self, n, k)
        self.plotResult(ax)
        dx  = (ax.get_xlim()[1] - ax.get_xlim()[0])/200.0
        dy  = (ax.get_ylim()[1] - ax.get_ylim()[0])/200.0
        self.x = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], dx)
        self.y = np.arange(ax.get_ylim()[0], ax.get_ylim()[1], dy)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.cov = np.array([[0.5,0.2],[0.2,0.5]])
        mu  = np.array(dpm_means())
        self.Z = biNormalDensity(mu[0], self.cov)(self.X, self.Y)
        for m in mu:
            self.Z = (self.Z + biNormalDensity(m, self.cov)(self.X,self.Y))
        manager = get_current_fig_manager()
        def updatefig(*args):
            try:
                ax.cla()
                self.sampleInteractively(1, ax)
                self.plotResult(ax)
                self.plotJoint(ax)
                manager.canvas.draw()
                return True
            except StopIteration:
                return False
        gobject.idle_add(updatefig)

    def sampleInteractively(self, n, ax):
        self.sample(n)
        mu  = np.array(dpm_means())
        for m in mu:
            self.Z = (self.Z + biNormalDensity(m, self.cov)(self.X,self.Y))
        self.Z/=2.0

def main():
    fig = figure()
    ax1 = fig.add_subplot(2,1,1, title="Data")
    ax2 = fig.add_subplot(2,1,2)
    dpm = InteractiveGDPM(40, 10, ax2)
    dpm.plotData(ax1)
    dpm.plotResult(ax2)
    show()

    fig = figure()
    ax  = fig.add_subplot(1,1,1, title="Statistics")
    dpm.plotStatistics(ax)
    show()

if __name__ == "__main__":
    sys.exit(main())
