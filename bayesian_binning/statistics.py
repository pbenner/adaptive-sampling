# Copyright (C) 2010 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import math
import numpy as np

def fac(n):
    return n-1 + abs(n-1) and fac(n-1)*n or 1L

def binomial(n,k):
    return fac(n) / (fac(k)*fac(n-k))

def binomialTransform(moments, n):
    result = math.pow(-moments[0], n)
    for k in range(1, n+1):
        result += binomial(n, k)*moments[k-1]*math.pow(-moments[0], n-k)
    return result

def centralMoments(moments_list, n):
    return map(lambda moments: binomialTransform(moments, n), zip(*moments_list))

def standardizedMoments(moments, n):
    m2 = centralMoments(moments, 2)
    mn = centralMoments(moments, n)
    return [ mu / math.pow(math.sqrt(var), n) for var, mu in zip(m2, mn) ]

def countStatistic(events):
    K = len(events)
    L = len(events[0])
    counts = np.zeros([K, L, L])
    for ks in range(0, L):
        for ke in range(ks, L):
            c = np.zeros(K)
            for i in range(ks, ke+1):
                for j in range(0, K):
                    c[j] += events[j][i]
            for j in range(0, K):
                counts[j][ks][ke] = c[j]
    return counts
