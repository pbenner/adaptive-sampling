# Copyright (C) 2011, 2012 Tobias Elze, Philipp Benner
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

#' Returns a default set of parameters for the binning and sampling algorithm.
#'
#' @param n.moments number of raw moments
#' @param model.posterior whether or not to compute the model posterior
#' @param bprob whether or not to compute the break probabilities
#' @param kl.psi compute the psi Kullback-Leibler divergence
#' @param kl.multibin compute the multibin Kullback-Leibler divergence
#' @param effective.counts compute the effective counts
#' @param effective.posterior.counts compute the effective posterior counts
#' @param density whether or not to compute the density
#' @param density.step step size for computing the density
#' @param density.range limits the range within which the density is computed
#' @param epsilon precision parameter for the extended prombs
#' @param threads number of threads that are used for computation
#' @param stacksize stacksize limit for multiple pthreads
#' @param algorithm 0: prombs, 1: multibin sampler
#' @param which specify the response for which all quantities are computed
#' @param hmm if 1 then hidden Markov models are used instead
#' @param rho cohesion parameter for the hidden Markov model 
#' @param samples the number of multibin samples for algorithm=1,
#' the first component of the vector specifies the number of burn-in samples
#' @examples
#' options <- make.options(model.posterior=0)
#' ls.str(options)
#' @export

make.options <-
  function(n.moments=2,
           model.posterior=1,
           bprob = 1,
           kl.psi=1,
           kl.multibin=0,
           effective.counts = 0,
           effective.posterior.counts = 0,
           density = 1,
           density.step = 0.01,
           density.range = c(0, 1),
           epsilon = 0.00001,
           threads = 1,
           stacksize = 256*1024,
           algorithm = 0,
           which = 0,
           hmm = 0,
           rho = 0.4,
           samples = c(100, 2000))
{
  env <- environment()
  env$n.moments                  <- n.moments
  env$model.posterior            <- model.posterior
  env$bprob                      <- bprob
  env$kl.psi                     <- kl.psi
  env$kl.multibin                <- kl.multibin
  env$effective.counts           <- effective.counts
  env$effective.posterior.counts <- effective.posterior.counts
  env$density                    <- density
  env$density.step               <- density.step
  env$density.range              <- density.range
  env$epsilon                    <- epsilon
  env$threads                    <- threads
  env$stacksize                  <- stacksize
  env$algorithm                  <- algorithm
  env$which                      <- which
  env$hmm                        <- hmm
  env$rho                        <- rho
  env$samples                    <- samples

  env
}
