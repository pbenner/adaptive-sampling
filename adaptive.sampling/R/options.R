# Copyright (C) 2011, 2012 Tobias Elze
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
