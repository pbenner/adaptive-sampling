#! /usr/bin/env python

import interface
import ConfigParser
import numpy as np
from matplotlib import *
from matplotlib.pyplot import *

config = ConfigParser.RawConfigParser()
config.read('../data/data2.cfg')

trials     = config.getint('Data', 'trials')
counts_str = config.get   ('Data', 'counts')
counts     = []

for value in counts_str.split(' '):
    counts.append(int(value))

result = interface.binning(counts,trials)

font = {'family'     : 'serif',
        'color'      : 'k',
        'weight' : 'normal',
        'size'   : 12,
        }

# x = np.arange(0, len(counts), 1)
# step(x, [a + b for a, b in zip(result[0], result[1])], 'r',
#      where='post', label='post')
# step(x, [a - b for a, b in zip(result[0], result[1])], 'r',
#      where='post', label='post')
# step(x, result[0], 'k', where='post', label='post')
# plot(x, [ float(a) / sum(counts) for a in counts], 'bo')
# xlabel('bin', font)
# ylabel('P(x)', font)
# show()

x = np.arange(0, len(counts), 1)
plot(x, [a + b for a, b in zip(result[0], result[1])], 'k--')
plot(x, [a - b for a, b in zip(result[0], result[1])], 'k--')
plot(x, result[0], 'r')
xlabel('bin', font)
ylabel('P(x)', font)
show()
