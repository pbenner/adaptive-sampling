#! /usr/bin/env python

# imports
################################################################################

import math

import numpy               as     np
import numpy.random.mtrand as     mt
import numpy.random        as     rd

from   matplotlib          import *
from   matplotlib.pyplot   import *
from   matplotlib.mlab     import *
from   matplotlib.image    import NonUniformImage
import matplotlib.patches  as     patches
import matplotlib.path     as     path

# statistics
################################################################################

def normalDensity(mu, sig2):
    """Returns the normal density"""
    stddev = np.sqrt(sig2)
    return (lambda x: mlab.normpdf(x, mu, stddev))

def biNormalDensity(mu, cov):
    return (lambda X, Y: mlab.bivariate_normal(X, Y, np.sqrt(cov[0,0]), np.sqrt(cov[1,1]), mu[0], mu[1], cov[0,1]))

def mNormalDensity(mu, cov):
    return (lambda x:
                np.power(np.linalg.det(cov),-0.5) *
            np.exp(-0.5*mu.shape[0]*np.log(2.0*np.pi)) *
            np.exp(-0.5*np.dot(np.dot(x-mu.transpose(),np.linalg.inv(cov)),x-mu)))

def randomElement(list, p):
    """Pick a random element from list with probabilities p"""
    val = rd.random()
    cp  = np.cumsum([0.0]+p)
    for i in range(0, len(list)):
        if cp[i] < val and val <= cp[i+1]:
            return list[i]
    return list[-1]
