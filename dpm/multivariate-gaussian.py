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

# local imports
################################################################################

from dpm        import *
from sampling   import *
from statistics import *

# Gaussian DPM
################################################################################

class MGaussianData(Data):
    def __init__(self):
        K    = 10
        N    = 20
        mu   = 50*rd.rand(K, 2)
        cov  = [[0.5,0.2],[0.2,0.5]]
        labeled_x = []
        for k in range(0, K):
            samples = rd.multivariate_normal(mu[k], cov, N)
            labeled_x.extend([ (elem,k) for elem in samples ])
        Data.__init__(self, labeled_x)
    def plot(self, ax):
        x = [ [] for i in range(0, self.N) ]
        for (x_,l_) in zip(self.x, self.labels):
            x[l_].append(x_)
        x = filter(lambda x_: not x_ == [], x)
        for x_ in x:
            ax.scatter(*zip(*x_), c=tuple(rd.rand(3)))

class MGaussianDPM(DPM):
    # parameters for the likelihood
    sig2   = [[0.5,0.2],[0.2,0.5]]
    # parameters for the prior
    mu_0   =  [1.0,1.0]
    sig2_0 = [[10.0,5.0],[5.0,10.0]]
    def __init__(self):
        data   = MGaussianData()
        DPM.__init__(self, data)
    def postPredDist(self, c):
        """Compute posterior predictive distribution"""
        data   = np.array(self.cl.getXByClass(c))
        sig2   = np.array(self.sig2)
        sig2_0 = np.array(self.sig2_0)
        sig2_inv   = np.linalg.inv(sig2)
        sig2_0_inv = np.linalg.inv(sig2_0)
        mu_0   = np.array(self.mu_0)
        num    = float(len(data))
        mean   = sum(data)/num
        sig2_n = np.linalg.inv(num*sig2_inv + sig2_0_inv)
        mu_n   = np.dot(sig2_n, np.dot(mu_0,sig2_0_inv) + num*np.dot(mean,sig2_inv))
        return mNormalDensity(mu_n, sig2_n + sig2)
    def predDist(self):
        """Compute predictive distribution"""
        sig2   = np.array(self.sig2)
        sig2_0 = np.array(self.sig2_0)
        mu_0   = np.array(self.mu_0)
        return mNormalDensity(mu_0, sig2_0 + sig2)
    def meanLikelihood(self):
        classes = self.cl.used_classes
        result  = 0
        for c in classes:
            x  = self.cl.getXByClass(c)
            mu = sum(x)/float(len(x))
            result += sum([ mNormalDensity(mu, self.sig2)(x_) for x_ in x ])/len(x)
        return result/len(classes)
    def plot(self, ax):
        x = [ self.cl.getXByClass(c) for c in self.cl.used_classes ]
        for x_ in x:
            ax.scatter(*zip(*x_), c=tuple(rd.rand(3)))
    def plotPrior(self, ax):
        dx, dy = 0.05, 0.05
        x = np.arange(0.0, 10.0, dx)
        y = np.arange(0.0, 10.0, dy)
        X,Y = np.meshgrid(x, y)
        Z = biNormalDensity(np.array(self.mu_0), np.array(self.sig2_0))(X, Y)
        im = NonUniformImage(ax, interpolation='bilinear', cmap=cm.gray)
        im.set_data(x, y, Z)

class MGaussianGibbsSampler(GibbsSampler):
    def printResult(self):
        print self.dpm.state()
        print self.dpm.da.labels
    def plotResult(self):
        fig = figure()
        ax1 = fig.add_subplot(2,1,1, title="Data")
        ax2 = fig.add_subplot(2,1,2, title="Clustering Result K="+str(len(self.dpm.cl.used_classes)))
        self.dpm.da.plot(ax1)
        self.dpm.plotPrior(ax2)
        self.dpm.plot(ax2)
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

# main
################################################################################

def main():
    dpm   = MGaussianDPM()
    gibbs = MGaussianGibbsSampler(dpm)
    gibbs.run(5)
    gibbs.plotResult()

if __name__ == "__main__":
    sys.exit(main())
