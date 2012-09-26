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

#' Calculate binning posterior quantities.
#' 
#' @param counts matrix of counts; each line one response option
#'  dimensions: K rows, L columns
#' @param alpha "pseudo counts"
#' @param beta relative class weights
#' @param gamma a priori importance of each consecutive bin
#' @param ... further options; see \code{\link{make.options}}
#' @references
#'  Poppe, S, Benner, P, Elze, T. 
#'  A predictive approach to nonparametric inference for adaptive 
#'  sequential sampling of psychophysical experiments.
#'  Journal of Mathematical Psychology 56 (2012) 179-195
#' @examples
#' L = 6 # number of stimuli
#' K = 2 # number of responses
#' counts.success <- c(2,3,2,4,7,7)
#' counts.failure <- c(8,7,7,6,3,2)
#' counts <- count.statistic(t(matrix(c(counts.success, counts.failure), L)))
#' alpha.success  <- c(1,1,1,1,1,1)
#' alpha.failure  <- c(1,1,1,1,1,1)
#' alpha  <- default.alpha(t(matrix(c(alpha.success, alpha.failure), L)))
#' beta   <- default.beta(L)
#' gamma  <- default.gamma(L)
#' result <- binning.posterior(counts, alpha, beta, gamma)
#' result <- binning.posterior(counts, alpha, beta, gamma, n.moments=5)
#' @export

binning.posterior <- function(counts, alpha, beta, gamma, ...) {
  L <- dim(counts)[1]
  K <- dim(counts)[3]
  storage.mode(counts) <- "double"
  storage.mode(alpha)  <- "double"
  storage.mode(beta)   <- "double"
  storage.mode(gamma)  <- "double"

  options <- make.options(...)

  marginal <- as.list(.Call("call_posterior",
                            counts, alpha, beta, gamma, options))
  attr(marginal, 'class') <- 'binning.posterior'
  marginal
}
