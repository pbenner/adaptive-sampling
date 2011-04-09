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
        self.dpm.history = []
        self.dpm.t       = []
    def run(self, n):
        """Run the Gibbs sampler for n iterations"""
        self.switches        = [0]
        self.mean_likelihood = [ ]
        self.mean_likelihood.append(self.dpm.meanLikelihood())
        self.dpm.history.append(self.dpm.cl)
        self.dpm.t.append(0)
        items = self.dpm.getItems()
        for i in range(0, n):
            ret = 0
            for item in items:
                ret += self.dpm.sample(item)
            self.switches.append(float(ret)/len(items))
            self.mean_likelihood.append(self.dpm.meanLikelihood())
            self.dpm.history.append(self.dpm.cl)
            self.dpm.t.append(self.dpm.t[-1]+1)
