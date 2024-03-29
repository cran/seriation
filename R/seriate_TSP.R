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

#' @import "TSP"

.tsp_control <- structure(
  list(
    method = "arbitrary insertion",
    rep = 10,
    two_opt = TRUE
  ),
  help = list(
    method = "used TSP method (see ? solve_TSP)",
    rep = "number of random restarts",
    two_opt = "use the 2-opt improvement heuristic?"
  )
)

seriate_dist_tsp <- function(x, control = NULL) {
  ## add a dummy city for cutting
  tsp <- insert_dummy(TSP(x), n = 1, label = "cut_here")

  if (is.null(control))
    control <- .tsp_control

  tour <- solve_TSP(tsp, method = control$method,
    control = control)

  o <- cut_tour(tour, cut = "cut_here", exclude_cut = TRUE)
  o
}

set_seriation_method(
  "dist",
  "TSP",
  seriate_dist_tsp,
  "Minimize Hamiltonian path length with a TSP solver.",
  .tsp_control,
  randomized = TRUE,
  optimizes = "Path_length"
)
