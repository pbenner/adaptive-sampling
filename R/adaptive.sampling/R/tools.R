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

#' Computes the full set of alpha parameters from a reduced set of pseudo counts.
#' 
#' @param alpha KxL matrix of pseudo counts where K is the number of
#' responses and L the number of stimuli
#' @examples
#' L = 6
#' alpha.success  <- c(1,1,1,1,1,1)
#' alpha.failure  <- c(1,1,1,1,1,1)
#' alpha  <- default.alpha(t(matrix(c(alpha.success, alpha.failure), L)))
#' @export

default.alpha <- function(alpha) {
  K <- dim(alpha)[1]
  L <- dim(alpha)[2]
  result <- array(dim=c(L,L,K))
  for (k in 1:K) {
    result[,,k] <- generate.alpha(alpha[k,])
  }
  result
}

generate.alpha <- function(alpha) {
  result <- outer(rep(1,length(alpha)),alpha)
  result[lower.tri(result)] <- 0
  for (i in 1:length(alpha)) {
    result[i,] <- cumsum(result[i,])
    for (j in i:length(alpha)) {
      result[i,j] <- result[i,j]/(j-i+1)
    }
  }
  result
}

#' Computes a set of default values for the beta parameter.
#' 
#' @param n number of stimuli
#' @export

default.beta <- function(n) {
  return(rep(1, n))
}

#' Computes a set of default values for the gamma parameter.
#' 
#' @param n number of stimuli
#' @export

default.gamma <- function(n) {
  result <- outer(rep(1,n),rep(1,n))
  result[lower.tri(result)] <- 0
  result
}

#' Computes the full count statistics from a set of bare counts.
#' 
#' @param counts KxL matrix of counts where K is the number of
#' responses and L the number of stimuli
#' @examples
#' L = 6 # number of stimuli
#' counts_success <- c(2,3,2,4,7,7)
#' counts_failure <- c(8,7,7,6,3,2)
#' counts <- count.statistic(t(matrix(c(counts_success, counts_failure), L)))
#' @export

count.statistic <- function(counts) {
  K <- dim(counts)[1]
  L <- dim(counts)[2]
  result <- array(dim=c(L,L,K))
  for (k in 1:K) {
    result[,,k] <- generate.statistic(counts[k,])
  }
  result
}

generate.statistic <- function(counts) {
  result <- outer(rep(1,length(counts)),counts)
  result[lower.tri(result)] <- 0
  for (i in 1:length(counts)) {
    result[i,] <- cumsum(result[i,])
  }
  result
}
