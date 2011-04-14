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
    def plotData(self, ax):
        num_clusters = dpm_num_clusters()
        for c in range(0, num_clusters):
            x, y = zip(*dpm_cluster(c))
            original_tags = dpm_original_tags(c)
            for (x_, y_, c_) in zip(x, y, original_tags):
                ax.scatter(x_, y_, c=self.cluster_colors[c_])
    def plotResult(self, ax):
        num_clusters = dpm_num_clusters()
#        print
#        print "number of clusters: "+str(num_clusters)
        for c in range(0, num_clusters):
#            print "plotting cluster: "+str(c)
#            print "cluster: "
#            print dpm_cluster(c)
#            if (len(dpm_cluster(c)) == 0):
#                dpm_print()
            x, y = zip(*dpm_cluster(c))
            ax.scatter(x, y, c=self.cluster_colors[c])
            ax.set_title("Clustering Result K="+str(num_clusters))

def main():
    dpm = GaussianDPM(10, 4)
    num_clusters = dpm.num_clusters()

    fig = figure()
    ax1 = fig.add_subplot(2,1,1, title="Data")
    ax2 = fig.add_subplot(2,1,2, title="Clustering Result K="+str(num_clusters))
    dpm.plotData(ax1)
    dpm.plotResult(ax2)

    manager = get_current_fig_manager()
    def updatefig(*args):
        try:
            dpm.sample(1)
            dpm.plotResult(ax2)
            manager.canvas.draw()
            return True
        except StopIteration:
            return False
    gobject.idle_add(updatefig)

    show()

if __name__ == "__main__":
    sys.exit(main())
