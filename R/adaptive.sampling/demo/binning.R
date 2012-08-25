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

require(adaptive.sampling)

binning.demo <- function() {
  L = 6 # number of stimuli
  K = 2 # number of responses
  counts_success <- c(2,3,2,4,7,7)
  counts_failure <- c(8,7,7,6,3,2)
  counts <- count.statistic(t(matrix(c(counts_success, counts_failure), L)))
  alpha_success  <- c(1,1,1,1,1,1)
  alpha_failure  <- c(1,1,1,1,1,1)
  alpha  <- default.alpha(t(matrix(c(alpha_success, alpha_failure), L)))
  beta   <- default.beta(L)
  gamma  <- default.gamma(L)

  binning.posterior(counts, alpha, beta, gamma)
}

marginal <- binning.demo()

plot(marginal)
