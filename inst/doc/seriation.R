## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)
library(seriation)
set.seed(1234)

## ----install, eval=FALSE------------------------------------------------------
# install.packages("seriation")

## ----load-package-------------------------------------------------------------
library(seriation)

## ----prepare-distance---------------------------------------------------------
data("SupremeCourt")
d <- as.dist(SupremeCourt)
d

## ----original-distance, fig.width=5, fig.height=5-----------------------------
pimage(d, main = "Original alphabetical order")

## ----find-order---------------------------------------------------------------
o <- seriate(d, method = "Spectral")
o

## ----inspect-order------------------------------------------------------------
get_order(o)

## ----apply-order--------------------------------------------------------------
d_ordered <- permute(d, o)
as.matrix(d_ordered)[1:4, 1:4]

## ----reordered-distance, fig.width=5, fig.height=5----------------------------
pimage(d, order = o, main = "Spectral seriation")

## ----assess-order-------------------------------------------------------------
rbind(
  original = criterion(d, method = c("2SUM", "Path_length")),
  seriated = criterion(d, o, method = c("2SUM", "Path_length"))
)

## ----matrix-order-------------------------------------------------------------
data("Wood")
dim(Wood)

o_matrix <- seriate(Wood, method = "Heatmap")
o_matrix

## ----matrix-permutations------------------------------------------------------
head(get_order(o_matrix, dim = 1))
get_order(o_matrix, dim = 2)

## ----matrix-images, fig.show="hold", out.width="49%", fig.width=5, fig.height=5----
pimage(Wood, main = "Original order")
pimage(Wood, order = o_matrix, main = "Seriated rows and columns")

## ----one-margin---------------------------------------------------------------
o_rows <- seriate(Wood, method = "PCA", margin = 1)
head(get_order(o_rows, dim = 1))

## ----methods------------------------------------------------------------------
head(list_seriation_methods("dist"))
list_seriation_methods("matrix")
get_seriation_method("dist", "Spectral")

