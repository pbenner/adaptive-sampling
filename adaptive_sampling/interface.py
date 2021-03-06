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

# allow to send ctrl-c to the library
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

# load interface
# ------------------------------------------------------------------------------

_lib = None

if   os.path.exists(os.path.dirname(__file__)+'/../src/.libs/libadaptive-sampling.so'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/.libs/libadaptive-sampling.so')
elif os.path.exists(os.path.dirname(__file__)+'/../src/.libs/libadaptive-sampling.dylib'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/.libs/libadaptive-sampling.dylib')
elif os.path.exists(os.path.dirname(__file__)+'/../src/.libs/cygadaptive-sampling-0.dll'):
     _lib = cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/.libs/cygadaptive-sampling-0.dll')
else:
     for libname in ['libadaptive-sampling.so.0', 'cygadaptive-sampling-0.dll', 'libadaptive-sampling.0.dylib']:
          if not _lib:
               try:
                    _lib = cdll.LoadLibrary(libname)
               except: pass

if not _lib:
     raise OSError('Couldn\'t find adaptive-sampling library.')

# structures
# ------------------------------------------------------------------------------

class VECTOR(Structure):
     _fields_ = [("size",     c_int),
                 ("content",  POINTER(c_double))]

class MATRIX(Structure):
     _fields_ = [("rows",    c_int),
                 ("columns", c_int),
                 ("content", POINTER(POINTER(c_double)))]

class DENSITY_RANGE(Structure):
     _fields_ = [("from", c_float),
                 ("to",   c_float)]

class OPTIONS(Structure):
     _fields_ = [("epsilon",              c_float),
                 ("verbose",              c_int),
                 ("prombsTest",           c_int),
                 ("bprob",                c_int),
                 ("threads",              c_int),
                 ("stacksize",            c_int),
                 ("kl_psi",               c_int),
                 ("kl_multibin",          c_int),
                 ("effective_counts",     c_int),
                 ("effective_posterior_counts", c_int),
                 ("which",                c_int),
                 ("algorithm",            c_int),
                 ("samples",            2*c_int),
                 ("density",              c_int),
                 ("density_step",         c_float),
                 ("density_range",        DENSITY_RANGE),
                 ("n_moments",            c_int),
                 ("n_density",            c_int),
                 ("model_posterior",      c_int),
                 ("hmm",                  c_int),
                 ("rho",                  c_float)]
     def __init__(self, options):
          self.which                = c_int(options["which"])
          self.threads              = c_int(options["threads"])
          self.stacksize            = c_int(options["stacksize"])
          self.samples[0]           = c_int(options["mgs_samples"][0])
          self.samples[1]           = c_int(options["mgs_samples"][1])
          self.density              = c_int(options["density"])
          self.density_step         = c_float(options["density_step"])
          self.density_range        = DENSITY_RANGE(*options["density_range"])
          self.n_moments            = c_int(options["n_moments"])
          self.n_density            = c_int(int(math.floor(1.0/options["density_step"]) + 1))
          self.epsilon              = c_float(options["epsilon"])
          self.verbose              = c_int(1) if options["verbose"]     else c_int(0)
          self.prombsTest           = c_int(1) if options["prombsTest"]  else c_int(0)
          self.bprob                = c_int(1) if options["bprob"]       else c_int(0)
          self.utility              = c_int(1) if options["utility"]     else c_int(0)
          self.kl_psi               = c_int(1) if options["kl_psi"]      else c_int(0)
          self.kl_multibin          = c_int(1) if options["kl_multibin"] else c_int(0)
          self.effective_counts     = c_int(1) if options["effective_counts"]  else c_int(0)
          self.effective_posterior_counts = c_int(1) if options["effective_posterior_counts"]  else c_int(0)
          self.model_posterior      = c_int(1) if options["model_posterior"]   else c_int(0)
          self.hmm                  = c_int(1) if options["hmm"]   else c_int(0)
          self.rho                  = c_float(options["rho"])
          if options["algorithm"] == "prombs":
               self.algorithm = c_int(0)
          elif options["algorithm"] == "mgs":
               self.algorithm = c_int(1)
          else:
               raise IOError("Unknown algorithm.")

class POSTERIOR(Structure):
     _fields_ = [("moments",   POINTER(MATRIX)),
                 ("density",   POINTER(MATRIX)),
                 ("bprob",     POINTER(VECTOR)),
                 ("mpost",     POINTER(VECTOR)),
                 ("utility",   POINTER(VECTOR))]

class UTILITY(Structure):
     _fields_ = [("expectation", POINTER(MATRIX)),
                 ("utility",     POINTER(VECTOR))]

# function prototypes
# ------------------------------------------------------------------------------

_lib._alloc_vector.restype   = POINTER(VECTOR)
_lib._alloc_vector.argtypes  = [c_int]

_lib._alloc_matrix.restype   = POINTER(MATRIX)
_lib._alloc_matrix.argtypes  = [c_int, c_int]

_lib._free_vector.restype    = None
_lib._free_vector.argtypes   = [POINTER(VECTOR)]

_lib._free_matrix.restype    = None
_lib._free_matrix.argtypes   = [POINTER(MATRIX)]

_lib._free.restype           = None
_lib._free.argtypes          = [POINTER(None)]

_lib._init_.restype          = None
_lib._init_.argtype          = [c_double]

_lib._free_.restype          = None
_lib._free_.argtype          = []

_lib.get_huge_val.restype    = c_double
_lib.get_huge_val.argtype    = []

_lib.gsl_sf_choose.restype   = c_double
_lib.gsl_sf_choose.argtype   = [c_uint, c_uint]

_lib.gsl_sf_lnchoose.restype = c_double
_lib.gsl_sf_lnchoose.argtype = [c_uint, c_uint]

_lib.posterior.restype       = POINTER(POSTERIOR)
_lib.posterior.argtypes      = [c_int, POINTER(POINTER(MATRIX)), POINTER(POINTER(MATRIX)), POINTER(VECTOR), POINTER(MATRIX), POINTER(OPTIONS)]

_lib.utility.restype         = POINTER(UTILITY)
_lib.utility.argtypes        = [c_int, POINTER(POINTER(MATRIX)), POINTER(POINTER(MATRIX)), POINTER(VECTOR), POINTER(MATRIX), POINTER(OPTIONS)]

_lib.utilityAt.restype       = POINTER(VECTOR)
_lib.utilityAt.argtypes      = [c_int, c_int, POINTER(POINTER(MATRIX)), POINTER(POINTER(MATRIX)), POINTER(VECTOR), POINTER(MATRIX), POINTER(OPTIONS)]

_lib.distance.restype        = c_double
_lib.distance.argtypes       = [c_int, c_int, c_int, POINTER(POINTER(MATRIX)), POINTER(POINTER(MATRIX)), POINTER(VECTOR), POINTER(MATRIX), POINTER(OPTIONS)]

# convert datatypes
# ------------------------------------------------------------------------------

def copyVectorToC(v, c_v):
     for i in range(0, c_v.contents.size):
          c_v.contents.content[i] = v[i]

def copyMatrixToC(m, c_m):
     for i in range(0, c_m.contents.rows):
          for j in range(0, c_m.contents.columns):
               c_m.contents.content[i][j] = m[i][j]

def getVector(c_v):
     v = []
     for i in range(0, c_v.contents.size):
          v.append(c_v.contents.content[i])
     return v

def getMatrix(c_m):
     m = []
     for i in range(0, c_m.contents.rows):
          m.append([])
          for j in range(0, c_m.contents.columns):
               m[i].append(c_m.contents.content[i][j])
     return m

# convert datatypes
# ------------------------------------------------------------------------------

def get_huge_val():
     return _lib.get_huge_val()

def gsl_sf_choose(n, m):
     return _lib.gsl_sf_choose(n, m)

def gsl_sf_lnchoose(n, m):
     return _lib.gsl_sf_lnchoose(n, m)

# 
# ------------------------------------------------------------------------------

def init(epsilon):
     c_epsilon = c_double(epsilon)
     _lib._init_(c_epsilon)

def free():
     _lib._free_()

def posterior(events, counts, alpha, beta, gamma, options):
     c_events = c_int(events)
     c_counts = (events*POINTER(MATRIX))()
     c_alpha  = (events*POINTER(MATRIX))()
     for i in range(0, events):
          c_counts[i]  = _lib._alloc_matrix(len(counts[i]), len(counts[i][0]))
          copyMatrixToC(counts[i], c_counts[i])
          c_alpha[i]   = _lib._alloc_matrix(len(alpha[i]), len(alpha[i][0]))
          copyMatrixToC(alpha[i],  c_alpha[i])
     c_beta  = _lib._alloc_vector(len(beta))
     copyVectorToC(beta,  c_beta)
     c_gamma = _lib._alloc_matrix(len(gamma), len(gamma[0]))
     copyMatrixToC(gamma,  c_gamma)
     c_options = pointer(OPTIONS(options))

     tmp = _lib.posterior(c_events, c_counts, c_alpha, c_beta, c_gamma, c_options)

     for i in range(0, events):
          _lib._free_matrix(c_counts[i])
          _lib._free_matrix(c_alpha[i])
     _lib._free_vector(c_beta)
     _lib._free_matrix(c_gamma)

     result = \
         { 'moments'   : getMatrix(tmp.contents.moments)   if bool(tmp.contents.moments)   else [],
           'density'   : getMatrix(tmp.contents.density)   if bool(tmp.contents.density)   else [],
           'bprob'     : getVector(tmp.contents.bprob)     if bool(tmp.contents.bprob)     else [],
           'mpost'     : getVector(tmp.contents.mpost)     if bool(tmp.contents.mpost)     else [] }

     if bool(tmp.contents.moments):
          _lib._free_matrix(tmp.contents.moments)
     if bool(tmp.contents.density):
          _lib._free_matrix(tmp.contents.density)
     if bool(tmp.contents.bprob):
          _lib._free_vector(tmp.contents.bprob)
     if bool(tmp.contents.mpost):
          _lib._free_vector(tmp.contents.mpost)
     _lib._free(tmp)

     return result

def utility(events, counts, alpha, beta, gamma, options):
     c_events      = c_int(events)
     c_counts      = (events*POINTER(MATRIX))()
     c_alpha       = (events*POINTER(MATRIX))()
     for i in range(0, events):
          c_counts[i]  = _lib._alloc_matrix(len(counts[i]), len(counts[i][0]))
          copyMatrixToC(counts[i], c_counts[i])
          c_alpha[i]   = _lib._alloc_matrix(len(alpha[i]), len(alpha[i][0]))
          copyMatrixToC(alpha[i],  c_alpha[i])
     c_beta  = _lib._alloc_vector(len(beta))
     copyVectorToC(beta,  c_beta)
     c_gamma = _lib._alloc_matrix(len(gamma), len(gamma[0]))
     copyMatrixToC(gamma,  c_gamma)
     c_options = pointer(OPTIONS(options))

     tmp = _lib.utility(c_events, c_counts, c_alpha, c_beta, c_gamma, c_options)

     for i in range(0, events):
          _lib._free_matrix(c_counts[i])
          _lib._free_matrix(c_alpha[i])
     _lib._free_vector(c_beta)
     _lib._free_matrix(c_gamma)

     result = \
         { 'expectation' : getMatrix(tmp.contents.expectation) if bool(tmp.contents.expectation) else [],
           'utility'     : getVector(tmp.contents.utility)     if bool(tmp.contents.utility)     else [] }

     if bool(tmp.contents.expectation):
          _lib._free_matrix(tmp.contents.expectation)
     if bool(tmp.contents.utility):
          _lib._free_vector(tmp.contents.utility)
     _lib._free(tmp)

     return result

def utilityAt(i, events, counts, alpha, beta, gamma, options):
     c_i           = c_int(i)
     c_events      = c_int(events)
     c_counts      = (events*POINTER(MATRIX))()
     c_alpha       = (events*POINTER(MATRIX))()
     for i in range(0, events):
          c_counts[i]  = _lib._alloc_matrix(len(counts[i]), len(counts[i][0]))
          copyMatrixToC(counts[i], c_counts[i])
          c_alpha[i]   = _lib._alloc_matrix(len(alpha[i]), len(alpha[i][0]))
          copyMatrixToC(alpha[i],  c_alpha[i])
     c_beta  = _lib._alloc_vector(len(beta))
     copyVectorToC(beta,  c_beta)
     c_gamma = _lib._alloc_matrix(len(gamma), len(gamma[0]))
     copyMatrixToC(gamma,  c_gamma)
     c_options = pointer(OPTIONS(options))

     c_result = _lib.utilityAt(c_i, c_events, c_counts, c_alpha, c_beta, c_gamma, c_options)
     result   = getVector(c_result)

     for i in range(0, events):
          _lib._free_matrix(c_counts[i])
          _lib._free_matrix(c_alpha[i])
     _lib._free_vector(c_beta)
     _lib._free_matrix(c_gamma)
     _lib._free_vector(c_result)

     return (result[0:-1], result[-1])

def distance(x, y, events, counts, alpha, beta, gamma, options):
     c_x           = c_int(x)
     c_y           = c_int(y)
     c_events      = c_int(events)
     c_counts      = (events*POINTER(MATRIX))()
     c_alpha       = (events*POINTER(MATRIX))()
     for i in range(0, events):
          c_counts[i]  = _lib._alloc_matrix(len(counts[i]), len(counts[i][0]))
          copyMatrixToC(counts[i], c_counts[i])
          c_alpha[i]   = _lib._alloc_matrix(len(alpha[i]), len(alpha[i][0]))
          copyMatrixToC(alpha[i],  c_alpha[i])
     c_beta  = _lib._alloc_vector(len(beta))
     copyVectorToC(beta,  c_beta)
     c_gamma = _lib._alloc_matrix(len(gamma), len(gamma[0]))
     copyMatrixToC(gamma,  c_gamma)
     c_options = pointer(OPTIONS(options))

     c_result = _lib.distance(c_x, c_y, c_events, c_counts, c_alpha, c_beta, c_gamma, c_options)

     for i in range(0, events):
          _lib._free_matrix(c_counts[i])
          _lib._free_matrix(c_alpha[i])
     _lib._free_vector(c_beta)
     _lib._free_matrix(c_gamma)

     return c_result
