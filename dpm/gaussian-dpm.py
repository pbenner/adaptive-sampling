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

from interface  import *

# Gaussian DPM
################################################################################

class GaussianDPM():
    def __init__(self, n, k):
        dpm_init(n, k)
        self.cluster_colors = [ tuple(rd.rand(3)) for i in range(0, n*k) ]
    def print_clusters(self):
        dpm_print()
    def num_clusters(self):
        return dpm_num_clusters()
    def sample(self, n):
        dpm_sample(n)
    def plot(self, ax):
        num_clusters = dpm_num_clusters()
        for c in range(0, num_clusters):
            x, y = zip(*dpm_cluster(c))
            print dpm_original_tags(c)
            ax.scatter(x, y, c=self.cluster_colors[c])

def main():
    dpm = GaussianDPM(100, 4)
    num_clusters = dpm.num_clusters()

    fig = figure()
    ax1 = fig.add_subplot(2,1,1, title="Data")
    ax2 = fig.add_subplot(2,1,2, title="Clustering Result K="+str(num_clusters))
    dpm.plot(ax2)
    show()

if __name__ == "__main__":
    sys.exit(main())
