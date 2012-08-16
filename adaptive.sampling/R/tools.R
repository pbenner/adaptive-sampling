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

default.beta <- function(n) {
  return(rep(1, n))
}

default.gamma <- function(n) {
  result <- outer(rep(1,n),rep(1,n))
  result[lower.tri(result)] <- 0
  result
}

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
