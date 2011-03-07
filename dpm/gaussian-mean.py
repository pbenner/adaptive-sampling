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

class GaussianData(Data):
    def __init__(self):
        K    = 10
        sig2 = 0.2
        N    = 50
        al   = np.cumsum(rd.randint(2,6,10))
        mu   = 100*mt.dirichlet(al, size=1)[0]
        labeled_x = []
        for k in range(0, K):
            samples = mt.normal(loc=mu[k], scale=np.sqrt(sig2), size=N)
            labeled_x.extend([ (elem,k) for elem in samples ])
        Data.__init__(self, labeled_x)
    def plotHist(self, ax):
        x = [ [] for i in range(0, self.N) ]
        for (x_,l_) in zip(self.x, self.labels):
            x[l_].append(x_)
        x = filter(lambda x_: not x_ == [], x)
        ax.hist(x, bins=self.N/10, histtype='barstacked')

class GaussianDPM(DPM):
    # parameters for the likelihood
    sig2   = 0.2
    # parameters for the prior
    mu_0   = 5.0
    sig2_0 = 0.3
    def __init__(self):
        data   = GaussianData()
        DPM.__init__(self, data)
    def postPredDist(self, c):
        """Compute posterior predictive distribution"""
        data   = self.cl.getXByClass(c)
        sig2   = self.sig2
        sig2_0 = self.sig2_0
        mu_0   = self.mu_0
        num    = float(len(data))
        mean   = sum(data)/num
        sig2_n = 1.0/(num/sig2 + 1.0/sig2_0)
        mu_n   = sig2_n*(mu_0/sig2_0 + num*mean/sig2)
        return normalDensity(mu_n, sig2_n + sig2)
    def predDist(self):
        """Compute predictive distribution"""
        return normalDensity(self.mu_0, self.sig2_0 + self.sig2)
    def meanLikelihood(self):
        classes = self.cl.used_classes
        result  = 0
        for c in classes:
            x  = self.cl.getXByClass(c)
            mu = float(sum(x))/len(x)
            result += sum([ normalDensity(mu, self.sig2)(x_) for x_ in x ])/len(x)
        return result/len(classes)
    def plotHist(self, ax):
        x = [ self.cl.getXByClass(c) for c in self.cl.used_classes ]
        ax.hist(x, bins=self.da.N/10, histtype='barstacked')

class GaussianGibbsSampler(GibbsSampler):
    def printResult(self):
        print self.dpm.state()
        print self.dpm.da.labels
    def plotResult(self):
        fig = figure()
        ax1 = fig.add_subplot(2,1,1, title="Data")
        ax2 = fig.add_subplot(2,1,2, title="Clustering Result")
        self.dpm.da.plotHist(ax1)
        self.dpm.plotHist(ax2)
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
    dpm   = GaussianDPM()
    gibbs = GaussianGibbsSampler(dpm)
    gibbs.run(100)
    gibbs.plotResult()

if __name__ == "__main__":
    sys.exit(main())
