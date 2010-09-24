#! /usr/bin/env python

import interface
import ConfigParser

config = ConfigParser.RawConfigParser()
config.read('../data/data1.cfg')

trials = config.getint('Data', 'trials')
counts = config.get('Data', 'counts')
x      = []

for value in counts.split(' '):
    x.append(int(value))

print interface.test([[1,2],[3,4]])
print interface.binning(x,trials)
print trials
print x
