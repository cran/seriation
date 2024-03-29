#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
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
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#' @rdname VAT
#' @export
ggVAT <- function(x,
  upper_tri = TRUE,
  lower_tri = TRUE,
  ...) {
  if (!inherits(x, "dist"))
    stop("x needs to be of class 'dist'!")
  ggpimage(x,
    seriate(x, "VAT"),
    upper_tri = upper_tri,
    lower_tri = lower_tri,
    ...)
}

#' @rdname VAT
#' @export
ggiVAT <- function(x,
  upper_tri = TRUE,
  lower_tri = TRUE,
  ...) {
  if (!inherits(x, "dist"))
    stop("x needs to be of class 'dist'!")
  x <- path_dist(x)
  ggpimage(x,
    seriate(x, "VAT"),
    upper_tri = upper_tri,
    lower_tri = lower_tri,
    ...)
}
