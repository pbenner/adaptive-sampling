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
            samples = mt.normal(loc=mu[k], scale=sig2, size=N)
            labeled_x.extend([ (elem,k) for elem in samples ])
        Data.__init__(self, labeled_x)

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
    def sample(self, item):
        """Perform a single sampling step"""
        old_class = self.cl.getClass(item)
        self.cl.release(item)
        x_i     = self.da.x[item]
        pred    = self.predDist()
        weights = []
        for c in self.cl.used_classes:
            postpred = self.postPredDist(c)
            weight   = float(self.cl.numItems(c))*postpred(x_i)
            weights.append(weight)
        weights.append(float(self.alpha)*pred(x_i))
        weights_norm = sum(weights)
        weights      = map(lambda w: w/weights_norm, weights)
        classes      = self.cl.used_classes+[self.cl.free_classes[0]]
        new_class    = randomElement(classes, weights)
        self.cl.assign(item, new_class)
        return not new_class == old_class

# main
################################################################################

def main():
    dpm   = GaussianDPM()
    gibbs = GibbsSampler(dpm)
    gibbs.run(200)
    gibbs.plotResult()

if __name__ == "__main__":
    sys.exit(main())
