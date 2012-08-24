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

require(prombs)

prombsExtended.demo <- function(m = 0) {
  g = c(1, 2, 3, 4, 5)
  f = c(1, 2, 3, 4, 5,
        0, 1, 2, 3, 4,
        0, 0, 1, 2, 3,
        0, 0, 0, 1, 2,
        0, 0, 0, 0, 1)
  h = c(1, 2, 3, 4, 5,
        0, 1, 2, 3, 4,
        0, 0, 1, 2, 3,
        0, 0, 0, 1, 2,
        0, 0, 0, 0, 1)

  dim(f) <- c(5,5)
  f      <- t(f)

  dim(h) <- c(5,5)
  h      <- t(h)

  # correct answer is [25.0063 200.0500 315.0788 160.0400 25.0063]
  #
  exp(prombsExtended(log(g), log(f), h, 0.0001, m));
}

result <- prombsExtended.demo()
