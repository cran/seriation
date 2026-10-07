## -----------------------------------------------------------------------------
library("seriation")
data("Wood")
Wood <- Wood[sample(nrow(Wood)), sample(ncol(Wood))]
dim(Wood)

DT::datatable(round(Wood, 2))

## -----------------------------------------------------------------------------
pimage(Wood)

## -----------------------------------------------------------------------------
o <- seriate(Wood, method = "Heatmap", seriation_method = "HC_Mean")
o

## -----------------------------------------------------------------------------
get_order(o, 2)

## -----------------------------------------------------------------------------
plot(o[[2]])
pimage(Wood, order = o)

## -----------------------------------------------------------------------------
pimage(Wood, order = "Heatmap", seriation_method = "HC_complete", 
       main = "Wood (hierarchical clustering)")
pimage(Wood, order = "Heatmap", seriation_method = "HC_Mean", 
       main = "Wood (reorder by row/col mean)")
pimage(Wood, order = "Heatmap", seriation_method = "GW_complete", 
       main = "Wood (reorder by Gruvaeus and Wainer heuristic)")
pimage(Wood, order = "Heatmap", 
       main = "Wood (default - optimal leaf ordering)")

## -----------------------------------------------------------------------------
hmap(Wood, method = "HC_complete", main = "Wood (hierarchical clustering)")
hmap(Wood, method = "HC_Mean", main = "Wood (reorder by row/col mean)")
hmap(Wood, method = "GW_complete", main = "Wood (reorder by Gruvaeus and Wainer heuristic)")
hmap(Wood, method = "OLO_complete", main = "Wood (opt. leaf ordering)")

## -----------------------------------------------------------------------------
register_DendSer()

hmap(Wood, method = "DendSer_BAR", main = "Wood (banded anti-Robinson)")

## -----------------------------------------------------------------------------
hmap(Wood, method = "HC_complete", 
     plot_margins = "distances",
     main = "Wood (hierarchical clustering)")

## -----------------------------------------------------------------------------
hmap(Wood, method = "MDS", main = "Wood (MDS)")
hmap(Wood, method = "MDS_angle", main = "Wood (Angle in 2D MDS space)")
hmap(Wood, method = "R2E", main = "Wood (Rank 2 ellipse seriation)")
hmap(Wood, method = "TSP", main = "Wood (Traveling salesperson)")

## -----------------------------------------------------------------------------
hmap(Wood, col = grays())
hmap(Wood, col = greenred())

hmap(Wood, col = colorRampPalette(c("brown", "orange", "red"))( 100 ) )

hmap(Wood, col = viridis::viridis(100))

## -----------------------------------------------------------------------------
if (!require("ggplot2")) install.packages("ggplot2")

library(ggplot2)
gghmap(Wood, method = "OLO")

## -----------------------------------------------------------------------------
o <- seriate(Wood, method = "Heatmap", seriation_method = "OLO")
heatmap(Wood, Rowv = as.dendrogram(o[[1]]), Colv = as.dendrogram(o[[2]]))

## -----------------------------------------------------------------------------
o <- seriate(Wood, method = "Heatmap", seriation_method = "Spectral")
heatmap(Wood, Rowv =  get_rank(o, 1), Colv =  get_rank(o, 2))

## -----------------------------------------------------------------------------
if (!suppressMessages(require("heatmaply"))) install.packages("heatmaply")

library("heatmaply")
heatmaply(Wood, seriate = "none", main = "HC")
heatmaply(Wood, seriate = "OLO", main = "OLO")

## -----------------------------------------------------------------------------
o <- seriate(Wood, method = "Heatmap", seriation_method = "OLO_ward")
heatmaply(Wood, Rowv = o[[1]], Colv = o[[2]], main = "OLO (Ward)")

o <- seriate(Wood, method = "Heatmap", seriation_method = "Spectral")
heatmaply(Wood, Rowv = get_rank(o, 1), Colv = get_rank(o, 2), main = "Spectral")

