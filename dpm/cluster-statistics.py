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

from dpm        import *
from sampling   import *
from statistics import *

# Gaussian DPM
################################################################################

class MGaussianData(Data):
    def __init__(self):
        K    = 2
        N    = 100
        mu   = [[1,1],[5,1]]
        cov  = [[0.5,0.2],[0.2,0.5]]
        labeled_x = []
        for k in range(0, K):
            samples = rd.multivariate_normal(mu[k], cov, N)
            labeled_x.extend([ (elem,k) for elem in samples ])
        Data.__init__(self, labeled_x)
        self.K   = K
        self.mu  = np.array(mu)
        self.cov = np.array(cov)
    def plot(self, ax):
        x = [ [] for i in range(0, self.N) ]
        for (x_,l_) in zip(self.x, self.labels):
            x[l_].append(x_)
        x = filter(lambda x_: not x_ == [], x)
        for x_ in x:
            xp, yp = zip(*x_)
            ax.scatter(xp, yp, c=tuple(rd.rand(3)))

        dx  = (ax.get_xlim()[1] - ax.get_xlim()[0])/200.0
        dy  = (ax.get_ylim()[1] - ax.get_ylim()[0])/200.0
        x   = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], dx)
        y   = np.arange(ax.get_ylim()[0], ax.get_ylim()[1], dy)
        X,Y = np.meshgrid(x, y)
        Z   = biNormalDensity(self.mu[0], self.cov)(X,Y)
        for k in range(1, self.K):
            Z = Z + biNormalDensity(self.mu[k], self.cov)(X,Y)
        im  = NonUniformImage(ax, interpolation='bilinear', cmap=cm.gray)
        im.set_data(x, y, Z)
        ax.images.append(im)

class MGaussianDPM(DPM):
    # parameters for the likelihood
    cov   = np.array([[0.5,0.2],[0.2,0.5]])
#    cov   = np.array([[0.1,0.01],[0.01,0.1]])
    # parameters for the prior
    mu_0  = np.array( [10.0,10.0])
#    cov_0 = np.array([[0.5,0.2],[0.2,0.5]])
    cov_0 = np.array([[10.0,5.0],[5.0,10.0]])
    def __init__(self):
        data   = MGaussianData()
        DPM.__init__(self, data)
    def postPredDist(self, c):
        """Compute posterior predictive distribution"""
        data   = np.array(self.cl.getXByClass(c))
        cov    = self.cov
        cov_0  = self.cov_0
        cov_inv   = np.linalg.inv(cov)
        cov_0_inv = np.linalg.inv(cov_0)
        mu_0   = self.mu_0
        num    = float(len(data))
        mean   = sum(data)/num
        cov_n  = np.linalg.inv(num*cov_inv + cov_0_inv)
        mu_n   = np.dot(cov_n, np.dot(mu_0,cov_0_inv) + num*np.dot(mean,cov_inv))
        return mNormalDensity(mu_n, cov_n + cov)
    def predDist(self):
        """Compute predictive distribution"""
        cov   = self.cov
        cov_0 = self.cov_0
        mu_0   = self.mu_0
        return mNormalDensity(mu_0, cov_0 + cov)
    def meanLikelihood(self):
        classes = self.cl.used_classes
        result  = 0
        for c in classes:
            x  = self.cl.getXByClass(c)
            mu = sum(x)/float(len(x))
            result += sum([ mNormalDensity(mu, self.cov)(x_) for x_ in x ])/len(x)
        return result/len(classes)
    def plot(self, ax):
        x = [ self.cl.getXByClass(c) for c in self.cl.used_classes ]
        for x_ in x:
            xp, yp = zip(*x_)
            ax.scatter(xp, yp, c=tuple(rd.rand(3)))
    def plotPrior(self, ax):
        dx  = (ax.get_xlim()[1] - ax.get_xlim()[0])/200.0
        dy  = (ax.get_ylim()[1] - ax.get_ylim()[0])/200.0
        x   = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], dx)
        y   = np.arange(ax.get_ylim()[0], ax.get_ylim()[1], dy)
        X,Y = np.meshgrid(x, y)
        Z   = biNormalDensity(self.mu_0, self.cov_0)(X,Y)
        im  = NonUniformImage(ax, interpolation='bilinear', cmap=cm.gray)
        im.set_data(x, y, Z)
        ax.images.append(im)

class MGaussianGibbsSampler(GibbsSampler):
    def printResult(self):
        print self.dpm.state()
        print self.dpm.da.labels
    def plotResult(self):
        fig = figure()
        ax1 = fig.add_subplot(2,1,1, title="Data")
        ax2 = fig.add_subplot(2,1,2, title="Clustering Result K="+str(len(self.dpm.cl.used_classes)))
        self.dpm.da.plot(ax1)
        self.dpm.cl.plot(ax2)
#        self.dpm.plotPrior(ax2)
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
    def plotAnimated(self):
        fig = figure()
        ax1 = fig.add_subplot(1,1,1, title="Statistics")
        ax2 = ax1.twinx()
        p1  = ax1.plot(self.switches)
        p2  = ax2.plot(self.mean_likelihood, color='green')
        ax1.set_ylabel("Mean class switches")
        ax2.set_ylabel("Mean likelihood")
        ax1.set_xlabel("iteration")
        ax1.legend([p1, p2], ["Mean class switches", "Mean likelihood"])
        fig = figure()
        ax1 = fig.add_subplot(2,1,1, title="Data")
        ax2 = fig.add_subplot(2,1,2, title="Clustering Result K="+str(len(self.dpm.cl.used_classes)))
        self.dpm.da.plot(ax1)
        self.dpm.cl.plot(ax2)
        manager = get_current_fig_manager()
        cli = iter(self.dpm.history)
        ti  = iter(self.dpm.t)
        def updatefig(*args):
            try:
                cli.next().plot(ax2)
                ax2.set_title("Clustering Result K="+str(len(self.dpm.cl.used_classes))+", t="+str(ti.next()))
                manager.canvas.draw()
                return True
            except StopIteration:
                return False
        gobject.idle_add(updatefig)
        show()

# main
################################################################################

def main():
    dpm   = MGaussianDPM()
    gibbs = MGaussianGibbsSampler(dpm)
    gibbs.run(100)
#    gibbs.plotResult()
    gibbs.plotAnimated()
 
if __name__ == "__main__":
    sys.exit(main())
