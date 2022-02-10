#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2015 Michael Hahsler, Christian Buchta and Kurt Hornik
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

## register seriation based on umap

register_umap <- function() {
  check_installed("umap")

  .contr <- unclass(umap::umap.defaults)
  .contr$n_components <- 1
  .contr$input <- "dist"

  umap_order <- function(x, control) {
    control <- .get_parameters(control, .contr)

    class(control) <- class(umap::umap.defaults)

    embedding <- umap::umap(as.matrix(x), config = control)
    order(embedding$layout)
  }

  set_seriation_method(
    "dist",
    "umap",
    umap_order,
    "Use 1D Uniform manifold approximation and projection (UMAP) embedding to create an order",
    .contr
  )
}