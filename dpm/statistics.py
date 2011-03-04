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

def randomElement(list, p):
    """Pick a random element from list with probabilities p"""
    val = rd.random()
    cp  = np.cumsum([0.0]+p)
    for i in range(0, len(list)):
        if cp[i] < val and val <= cp[i+1]:
            return list[i]
    return list[-1]
