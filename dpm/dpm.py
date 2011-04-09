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

from   statistics import *

# general data structures
################################################################################

class Data():
    def __init__(self, labeled_x):
        rd.shuffle(labeled_x)
        self.x      = [ elem  for (elem,index) in labeled_x ]
        self.labels = [ index for (elem,index) in labeled_x ]
        self.N      = len(self.x)
        self.items  = range(0, self.N)
    def __iter__(self):
        return iter(self.x)
    def __getitem__(self, index):
        return self.x[index]

class ClassAssignments():
    _INIT_NUM_CLASSES = 20
    def __init__(self, data):
        self.data = data
        self.classes           = [ []    for i in range(0, data.N) ]
        self.class_assignments = [ None  for i in range(0, data.N) ]
        self.used_classes      = []
        self.free_classes      = range(0, data.N)
        # randomly assign all items to some class
        for item in range(0, data.N):
            self.assign(item, rd.randint(0, self._INIT_NUM_CLASSES))
    def getXByClass(self, c):
        """Return all x values of a class"""
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
        if len(self.classes[new_class]) == 0:
            self.used_classes.append(new_class)
            self.free_classes.remove(new_class)
        self.classes[new_class].append(item)
        self.class_assignments[item] = new_class
    def reassign(self, item, new_class):
        """Reassign an item to some other class"""
        self.release(item)
        self.assign(item, new_class)
    def plot(self, ax):
        x = [ self.getXByClass(c) for c in self.used_classes ]
        for x_ in x:
            xp, yp = zip(*x_)
            ax.scatter(xp, yp, c=tuple(rd.rand(3)))

class DPM():
    # parameters for the dirichlet process
    alpha  = 1.0
    def __init__(self, data):
        self.da = data
        self.cl = ClassAssignments(data)
    def state(self):
        return self.cl.class_assignments
    def getItems(self):
        return self.da.items
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
