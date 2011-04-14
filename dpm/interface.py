# Copyright (C) 2011 Philipp Benner
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

import os
import numpy as np
import math

from ctypes import *

# load interface
# ------------------------------------------------------------------------------

_lib = None

if   os.path.exists(os.path.dirname(__file__)+'/.libs/libdpm.so'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libdpm.so')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/libdpm.dylib'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libdpm.dylib')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/cygdpm-0.dll'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/cygdpm-0.dll')
else:
     for libname in ['libdpm.so.0', 'cygdpm-0.dll', 'libdpm.0.dylib']:
          if not _lib:
               try:
                    _lib = cdll.LoadLibrary(libname)
               except: pass

if not _lib:
     raise OSError('Couldn\'t find dpm library.')

# structures
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size", c_int),
                 ("vec",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_int),
                 ("columns", c_int),
                 ("mat",     POINTER(POINTER(c_double)))]

class MARGINAL_RANGE(Structure):
     _fields_ = [("from", c_float),
                 ("to",   c_float)]

class OPTIONS(Structure):
     _fields_ = [("epsilon",           c_float),
                 ("verbose",           c_int),
                 ("prombsTest",        c_int),
                 ("bprob",             c_int),
                 ("differential_gain", c_int),
                 ("effective_counts",  c_int),
                 ("multibin_entropy",  c_int),
                 ("which",             c_int),
                 ("marginal",          c_int),
                 ("marginal_step",     c_float),
                 ("marginal_range",    MARGINAL_RANGE),
                 ("n_moments",         c_int),
                 ("n_marginals",       c_int),
                 ("model_posterior",   c_int)]
     def __init__(self, options):
          self.which          = c_int(options["which"])
          self.marginal       = c_int(options["marginal"])
          self.marginal_step  = c_float(options["marginal_step"])
          self.marginal_range = MARGINAL_RANGE(*options["marginal_range"])
          self.n_moments      = c_int(options["n_moments"])
          self.n_marginals    = c_int(int(math.floor(1.0/options["marginal_step"]) + 1))
          self.epsilon        = c_float(options["epsilon"])
          self.verbose        = c_int(1) if options["verbose"]    else c_int(0)
          self.prombsTest     = c_int(1) if options["prombsTest"] else c_int(0)
          self.bprob          = c_int(1) if options["bprob"]      else c_int(0)
          self.differential_gain = c_int(1) if options["differential_gain"] else c_int(0)
          self.effective_counts  = c_int(1) if options["effective_counts"]  else c_int(0)
          self.multibin_entropy  = c_int(1) if options["multibin_entropy"]  else c_int(0)
          self.model_posterior   = c_int(1) if options["model_posterior"]   else c_int(0)

class BINNING_RESULT(Structure):
     _fields_ = [("moments",           POINTER(MATRIX)),
                 ("marginals",         POINTER(MATRIX)),
                 ("bprob",             POINTER(VECTOR)),
                 ("mpost",             POINTER(VECTOR)),
                 ("differential_gain", POINTER(VECTOR)),
                 ("effective_counts",  POINTER(VECTOR)),
                 ("multibin_entropy",  POINTER(VECTOR))]

# function prototypes
# ------------------------------------------------------------------------------

_lib._allocVector.restype  = POINTER(VECTOR)
_lib._allocVector.argtypes = [c_int]

_lib._allocMatrix.restype  = POINTER(MATRIX)
_lib._allocMatrix.argtypes = [c_int, c_int]

_lib._freeVector.restype   = None
_lib._freeVector.argtypes  = [POINTER(VECTOR)]

_lib._freeMatrix.restype   = None
_lib._freeMatrix.argtypes  = [POINTER(MATRIX)]

_lib._free.restype         = None
_lib._free.argtypes        = [POINTER(None)]

_lib._dpm_init.restype     = None
_lib._dpm_init.argtypes    = [c_uint, c_uint]

_lib._dpm_num_clusters.restype  = c_uint
_lib._dpm_num_clusters.argtypes = []

_lib._dpm_cluster.restype  = POINTER(MATRIX)
_lib._dpm_cluster.argtypes = [c_uint]

_lib._dpm_print.restype    = None
_lib._dpm_print.argtypes   = []

_lib._dpm_sample.restype   = None
_lib._dpm_sample.argtypes  = []

_lib._dpm_free.restype     = None
_lib._dpm_free.argtypes    = []

# convert datatypes
# ------------------------------------------------------------------------------

def copyVectorToC(v, c_v):
     for i in range(0, c_v.contents.size):
          c_v.contents.vec[i] = v[i]

def copyMatrixToC(m, c_m):
     for i in range(0, c_m.contents.rows):
          for j in range(0, c_m.contents.columns):
               c_m.contents.mat[i][j] = m[i][j]

def getVector(c_v):
     v = []
     for i in range(0, c_v.contents.size):
          v.append(c_v.contents.vec[i])
     return v

def getMatrix(c_m):
     m = []
     for i in range(0, c_m.contents.rows):
          m.append([])
          for j in range(0, c_m.contents.columns):
               m[i].append(c_m.contents.mat[i][j])
     return m

# 
# ------------------------------------------------------------------------------

def dpm_init(n, k):
     _lib._dpm_init(n, k)

def dpm_print():
     _lib._dpm_print()

def dpm_num_clusters():
     return _lib._dpm_num_clusters()

def dpm_cluster(c):
     result  = _lib._dpm_cluster(c)
     cluster = getMatrix(result)
     _lib._freeMatrix(result)
     return cluster

def dpm_sample(n):
     _lib_dpm_sample(n)

def dpm_free():
     _lib_dpm_free()
