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

from ctypes      import *
from ctypes.util import find_library

# load interface
# ------------------------------------------------------------------------------

if   os.path.exists(os.path.dirname(__file__)+'/.libs/libbayesian-binning.so'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libbayesian-binning.so')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/libbayesian-binning.dylib'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/libbayesian-binning.dylib')
elif os.path.exists(os.path.dirname(__file__)+'/.libs/cygbayesian-binning-0.dll'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/.libs/cygbayesian-binning-0.dll')
elif find_library('bayesian-binning'):
     _lib = cdll.LoadLibrary(find_library('bayesian-binning'))
else:
     try:
          _lib = cdll.LoadLibrary('libbayesian-binning.so.0')
     except OSError:
          _lib = cdll.LoadLibrary('cygbayesian-binning-0.dll')

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

_lib.binning.restype       = POINTER(BINNING_RESULT)
_lib.binning.argtypes      = [POINTER(MATRIX), POINTER(MATRIX), POINTER(VECTOR), POINTER(OPTIONS)]

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

def binning(counts, alpha, prior, options):
     c_counts  = _lib._allocMatrix(len(counts), len(counts[0]))
     copyMatrixToC(counts, c_counts)
     c_alpha   = _lib._allocMatrix(len(alpha), len(alpha[0]))
     copyMatrixToC(alpha,  c_alpha)
     c_prior   = _lib._allocVector(len(prior))
     copyVectorToC(prior,  c_prior)
     c_options = pointer(OPTIONS(options))

     tmp = _lib.binning(c_counts, c_alpha, c_prior, c_options)

     _lib._freeMatrix(c_counts)
     _lib._freeMatrix(c_alpha)
     _lib._freeVector(c_prior)

     result = \
     { 'moments'   : getMatrix(tmp.contents.moments)   if tmp.contents.moments   else [],
       'marginals' : getMatrix(tmp.contents.marginals) if tmp.contents.marginals else [],
       'bprob'     : getVector(tmp.contents.bprob),
       'mpost'     : getVector(tmp.contents.mpost),
       'differential_gain' : getVector(tmp.contents.differential_gain),
       'effective_counts'  : getVector(tmp.contents.effective_counts),
       'multibin_entropy'  : getVector(tmp.contents.multibin_entropy) }

     if tmp.contents.moments:
          _lib._freeMatrix(tmp.contents.moments)
     if tmp.contents.marginals:
          _lib._freeMatrix(tmp.contents.marginals)
     _lib._freeVector(tmp.contents.bprob)
     _lib._freeVector(tmp.contents.mpost)
     _lib._freeVector(tmp.contents.differential_gain)
     _lib._freeVector(tmp.contents.effective_counts)
     _lib._freeVector(tmp.contents.multibin_entropy)
     _lib._free(tmp)

     return result
