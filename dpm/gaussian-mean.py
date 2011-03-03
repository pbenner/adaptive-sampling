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

# statistics
################################################################################

def normalDensity(mu, sig2):
    """Returns the normal density"""
    return (lambda x: 1/(np.sqrt(sig2) * np.sqrt(2 * np.pi)) *
            np.exp( - (x - mu)**2 / (2 * np.sqrt(sig2)**2)))

def randomElem(list, p):
    """Pick a random element from list with probabilities p"""
    val = rd.random()
    cp  = np.cumsum([0.0]+p)
    for i in range(0, len(list)-1):
        if cp[i] < val and val <= cp[i+1]:
            return list[i]

# data structures
################################################################################

class Data():
    def __init__(self, labeled_x):
        rd.shuffle(labeled_x)
        self.x      = [ elem  for (elem,index) in labeled_x ]
        self.labels = [ index for (elem,index) in labeled_x ]
        self.N      = len(self.x)
    def __iter__(self):
        return iter(self.x)
    def __getitem__(self, index):
        return self.x[index]
    def plotHist(self):
        count, bins, ignored = hist(self.x, bins=self.N/10)
        show()

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

class ClassAssignments():
    _INIT_NUM_CLASSES = 10
    def __init__(self, data):
        self.data = data
        self.classes           = [ [] for i in range(0, data.N) ]
        self.class_assignments = [ 0  for i in range(0, data.N) ]
        self.used_classes      = range(0, self._INIT_NUM_CLASSES)
        self.free_classes      = range(self._INIT_NUM_CLASSES, data.N)
        # randomly assign all items to some class
        for item in range(0, data.N):
            self.assign(item, rd.randint(0, self._INIT_NUM_CLASSES))
    def getItemsByClass(self, c):
        """Return all items of a class"""
        return [ self.data[index] for index in self.classes[c] ]
    def getClass(self, item):
        """Get the class of an item"""
        return self.class_assignments[item]
    def numItems(self, c):
        """Count the number of items of a class"""
        return len(self.classes[c])
    def release(self, item):
        """Remove an item from its class"""
        old_class = self.getClass(item)
        self.classes[old_class].remove(item)
        self.class_assignments[item] = None
        if len(self.classes[old_class]) == 0:
            self.used_classes.remove(old_class)
            self.free_classes.append(old_class)
    def assign(self, item, new_class):
        """Assign an item to some class"""
        self.classes[new_class].append(item)
        self.class_assignments[item] = new_class
    def reassign(self, item, new_class):
        """Reassign an item to some other class"""
        self.release(item)
        self.assign(item, new_class)

class GaussianDPM():
    # parameters for the dirichlet process
    alpha  = 1.0
    # parameters for the likelihood
    sig2   = 0.2
    # parameters for the prior
    mu_0   = 5.0
    sig2_0 = 0.3
    def __init__(self):
        self.da = GaussianData()
        self.cl = ClassAssignments(self.da)
    def predictivePosterior(self, c):
        data   = self.cl.getItemsByClass(c)
        sig2   = self.sig2
        sig2_0 = self.sig2_0
        mu_0   = self.mu_0
        num    = float(len(data))
        mean   = sum(data)/num
        sig2_n = 1.0/(num/sig2 + 1.0/sig2_0)
        mu_n   = sig2_n*(mu_0/sig2_0 + num*mean/sig2)
        return mu_n, sig2_n + sig2
    def predictive(self):
        return self.mu_0, self.sig2_0 + self.sig2
    def sample(self, item):
        self.cl.release(item)
        x_i     = self.da.x[item]
        pred    = self.predictive()
        weights = []
        for c in self.cl.used_classes:
            predpost = self.predictivePosterior(c)
            weight   = (float(self.cl.numItems(c)) *
                        normalDensity(*predpost)(x_i))
            weights.append(weight)
        weights.append(float(self.alpha) *
                       normalDensity(*pred)(x_i))
        weights = map(lambda w: w/sum(weights), weights)
        print weights
        print sum(weights)

# main
################################################################################

da = GaussianData()
cl = ClassAssignments(da)
print cl.class_assignments
print cl.classes

#dpm = GaussianDPM()
#dpm.sample(1)
