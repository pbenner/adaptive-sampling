#! /usr/bin/env python

import os.path
import numpy as np
import ConfigParser

from itertools import izip

try:
    from matplotlib import *
    from matplotlib.pyplot import *
    from matplotlib.image import NonUniformImage
    import matplotlib.patches as patches
    import matplotlib.path as path
except ImportError:
    print "Error: Couldn't load matplotlib."
    exit(1)
except RuntimeError:
    print "Error: No x11 available."
    exit(1)

# tools
# ------------------------------------------------------------------------------

def argmax(array):
    result = []
    for i in izip(array, xrange(len(array))):
        if i[0] == max(array):
            result.append(i[1])
    return result

# functions for reading config options
# ------------------------------------------------------------------------------

def readVector(config, section, option, converter):
    vector_str = config.get(section, option)
    vector     = map(converter, vector_str.split(' '))
    return vector

# parse config
# ------------------------------------------------------------------------------

def parseConfig(config_file):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_file)

    return readVector(config_parser, 'Sampling Result', 'utility', float)

# plot entropies
# ------------------------------------------------------------------------------

font = {'family'     : 'serif',
        'weight'     : 'normal',
        'size'       : 12 }

def plotMaxima(ax, utility):
    x = map(lambda x: x+1, argmax(utility))
    if (len(x) < len(utility)):
        y = [ (ax.get_ylim()[0] + ax.get_ylim()[1])/2.0 ] * len(x)
        ax.plot(x,y, "k+", mew=3, ms=12)

def plotEntropies(targets, filename, title):
    fig = figure()
    fig.subplots_adjust(hspace=0.0, wspace=0.4)
    ax1 = fig.add_subplot(4,1,1, title=title)
    ax2 = fig.add_subplot(4,1,2)
    ax3 = fig.add_subplot(4,1,3)
    ax4 = fig.add_subplot(4,1,4)
    x = range(1, 36)
    ax1.plot(x, targets[0]['utility'], label=targets[0]['name'], color=targets[0]['color'], lw=2)
    plotMaxima(ax1, targets[0]['utility'])
    ax2.plot(x, targets[1]['utility'], label=targets[1]['name'], color=targets[1]['color'], lw=2)
    plotMaxima(ax2, targets[1]['utility'])
    ax3.plot(x, targets[2]['utility'], label=targets[2]['name'], color=targets[2]['color'], lw=2)
    plotMaxima(ax3, targets[2]['utility'])
    ax4.plot(x, targets[3]['utility'], label=targets[3]['name'], color=targets[3]['color'], lw=2)
    plotMaxima(ax4, targets[3]['utility'])
    ax4.set_xlabel('x')
    ax1.set_ylabel(targets[0]['name'], font)
    ax2.set_ylabel(targets[1]['name'], font)
    ax3.set_ylabel(targets[2]['name'], font)
    ax4.set_ylabel(targets[3]['name'], font)
    ax1.set_xlim(1, 35)
    ax2.set_xlim(1, 35)
    ax3.set_xlim(1, 35)
    ax4.set_xlim(1, 35)
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xticks([])
    ax1.set_yticks([])
    ax2.set_yticks([])
    ax3.set_yticks([])
    ax4.set_yticks([])
    savefig(filename, bbox_inches='tight', pad_inches=0.2)
#    savefig(filename)
#    show()

# main
# ------------------------------------------------------------------------------

targets1 = [
    { 'file': 'example2-tmp/example2.result01', 'color': 'blue',    'name': r'$U^*_1$' },
    { 'file': 'example2-tmp/example2.result02', 'color': 'red',     'name': r'$U^*_2$' },
    { 'file': 'example2-tmp/example2.result03', 'color': '#ee8d18', 'name': r'$U^*_3$' },
    { 'file': 'example2-tmp/example2.result04', 'color': 'green',   'name': r'$U^*_4$' }
    ]
targets2 = [
    { 'file': 'example2-tmp/example2.result05', 'color': 'blue',    'name': r'$U^*_1$' },
    { 'file': 'example2-tmp/example2.result06', 'color': 'red',     'name': r'$U^*_2$' },
    { 'file': 'example2-tmp/example2.result07', 'color': '#ee8d18', 'name': r'$U^*_3$' },
    { 'file': 'example2-tmp/example2.result08', 'color': 'green',   'name': r'$U^*_4$' }
    ]
targets3 = [
    { 'file': 'example2-tmp/example2.result09', 'color': 'blue',    'name': r'$U^*_1$' },
    { 'file': 'example2-tmp/example2.result10', 'color': 'red',     'name': r'$U^*_2$' },
    { 'file': 'example2-tmp/example2.result11', 'color': '#ee8d18', 'name': r'$U^*_3$' },
    { 'file': 'example2-tmp/example2.result12', 'color': 'green',   'name': r'$U^*_4$' }
    ]
targets4 = [
    { 'file': 'example2-tmp/example2.result13', 'color': 'blue',    'name': r'$U^*_1$' },
    { 'file': 'example2-tmp/example2.result14', 'color': 'red',     'name': r'$U^*_2$' },
    { 'file': 'example2-tmp/example2.result15', 'color': '#ee8d18', 'name': r'$U^*_3$' },
    { 'file': 'example2-tmp/example2.result16', 'color': 'green',   'name': r'$U^*_4$' }
    ]
targets5 = [
    { 'file': 'example2-tmp/example2.result17', 'color': 'blue',    'name': r'$U^*_1$' },
    { 'file': 'example2-tmp/example2.result18', 'color': 'red',     'name': r'$U^*_2$' },
    { 'file': 'example2-tmp/example2.result19', 'color': '#ee8d18', 'name': r'$U^*_3$' },
    { 'file': 'example2-tmp/example2.result20', 'color': 'green',   'name': r'$U^*_4$' }
    ]

def main():
    for target in targets1:
        target['utility'] = parseConfig(target['file'])
    for target in targets2:
        target['utility'] = parseConfig(target['file'])
    for target in targets3:
        target['utility'] = parseConfig(target['file'])
    for target in targets4:
        target['utility'] = parseConfig(target['file'])
    for target in targets5:
        target['utility'] = parseConfig(target['file'])
    plotEntropies(targets1, 'example2-tmp/example2-1.pdf', "(a) 0 Samples")
    plotEntropies(targets2, 'example2-tmp/example2-2.pdf', "(a) 1 Sample")
    plotEntropies(targets3, 'example2-tmp/example2-3.pdf', "(a) 3 Samples")
    plotEntropies(targets4, 'example2-tmp/example2-4.pdf', "(b) 6 Samples")
    plotEntropies(targets5, 'example2-tmp/example2-5.pdf', "(c) 200 Samples")

if __name__ == "__main__":
    sys.exit(main())
