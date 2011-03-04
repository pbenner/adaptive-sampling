#! /usr/bin/env python

# imports
################################################################################

import math

import numpy               as     np
import numpy.random.mtrand as     mt
import numpy.random        as     rd

from   matplotlib          import *
from   matplotlib.pyplot   import *
from   matplotlib.image    import NonUniformImage
import matplotlib.patches  as     patches
import matplotlib.path     as     path

# sampling
################################################################################

class GibbsSampler():
    def __init__(self, dpm):
        self.dpm = dpm
    def run(self, n):
        """Run the Gibbs sampler for n iterations"""
        self.switches        = []
        self.mean_likelihood = []
        items = self.dpm.getItems()
        for i in range(0, n):
            ret = 0
            for item in items:
                ret += self.dpm.sample(item)
            self.switches.append(float(ret)/len(items))
            self.mean_likelihood.append(self.dpm.meanLikelihood())
    def printResult(self):
        print self.dpm.state()
        print self.dpm.da.labels
    def plotResult(self):
        fig = figure()
        ax1 = fig.add_subplot(2,1,1, title="Data")
        ax2 = fig.add_subplot(2,1,2, title="Clustering Result")
        self.dpm.da.plotHist(ax1)
        self.dpm.cl.plotHist(ax2)
        fig = figure()
        ax1 = fig.add_subplot(1,1,1, title="Statistics")
        ax2 = ax1.twinx()
        p1  = ax1.plot(self.switches)
        p2  = ax2.plot(self.mean_likelihood, color='green')
        ax1.set_ylabel("Mean class switches")
        ax2.set_ylabel("Mean likelihood")
        ax1.set_xlabel("iteration")
        ax1.legend([p1, p2], ["Mean class switches", "Mean likelihood"])
        show()
