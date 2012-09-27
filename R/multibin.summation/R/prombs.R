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

#' Proximal MultiBin Summation (ProMBS)
#' 
#' @param g vector of length L with values on log scale
#' @param f LxL upper triangular matrix with values on log scale
#' @param m integer m <= L that specifies the number of products
#' @seealso \code{\link{prombsExtended}}
#' @references
#'  Poppe, S, Benner, P, Elze, T. 
#'  A predictive approach to nonparametric inference for adaptive 
#'  sequential sampling of psychophysical experiments.
#'  Journal of Mathematical Psychology 56 (2012) 179-195
#' @examples
#' g = c(1, 2, 3, 4, 5)
#' f = c(1, 2, 3, 4, 5,
#'       0, 1, 2, 3, 4,
#'       0, 0, 1, 2, 3,
#'       0, 0, 0, 1, 2,
#'       0, 0, 0, 0, 1)
#'
#' dim(f) <- c(5,5)
#' f      <- t(f)
#'
#' exp(prombs(log(g), log(f)));
#' @export

prombs <- function(g, f, m=0) {
  L <- length(g)

  if (m == 0) {
    m <- L
  }
  dim(g) <- c(L)
  dim(f) <- c(L, L)
  dim(m) <- c(1)

  storage.mode(g) <- "double"
  storage.mode(f) <- "double"
  storage.mode(m) <- "integer"

  .Call("call_prombs", g, f, m, "prombs")
}

#' Extended Proximal MultiBin Summation (ProMBS)
#' 
#' @param g vector of length L with values on log scale
#' @param f LxL upper triangular matrix with values on log scale
#' @param h LxL upper triangular matrix with values on normal scale
#' @param epsilon precision parameter
#' @param m integer m <= L that specifies the number of products
#' @seealso \code{\link{prombs}}
#' @references
#'  Poppe, S, Benner, P, Elze, T. 
#'  A predictive approach to nonparametric inference for adaptive 
#'  sequential sampling of psychophysical experiments.
#'  Journal of Mathematical Psychology 56 (2012) 179-195
#' @examples
#' g = c(1, 2, 3, 4, 5)
#' f = c(1, 2, 3, 4, 5,
#'       0, 1, 2, 3, 4,
#'       0, 0, 1, 2, 3,
#'       0, 0, 0, 1, 2,
#'       0, 0, 0, 0, 1)
#' h = c(1, 2, 3, 4, 5,
#'       0, 1, 2, 3, 4,
#'       0, 0, 1, 2, 3,
#'       0, 0, 0, 1, 2,
#'       0, 0, 0, 0, 1)
#'
#' dim(f) <- c(5,5)
#' f      <- t(f)
#'
#' dim(h) <- c(5,5)
#' h      <- t(h)
#'
#' exp(prombsExtended(log(g), log(f), h, 0.0001));
#' @export

prombsExtended <- function(g, f, h, epsilon, m=0) {
  L <- length(g)

  if (m == 0) {
    m <- L
  }
  dim(g)       <- c(L)
  dim(f)       <- c(L, L)
  dim(h)       <- c(L, L)
  dim(epsilon) <- c(1)
  dim(m)       <- c(1)

  storage.mode(g)       <- "double"
  storage.mode(f)       <- "double"
  storage.mode(h)       <- "double"
  storage.mode(m)       <- "integer"
  storage.mode(epsilon) <- "double"

  .Call("call_prombs_extended", g, f, h, epsilon, m, "prombs")
}
