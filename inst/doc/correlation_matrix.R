## -----------------------------------------------------------------------------
if (!require("seriation")) install.packages("seriation")

library("seriation")
data("mtcars")

DT::datatable(mtcars)

## -----------------------------------------------------------------------------
m <- cor(mtcars)
round(m, 2)

## -----------------------------------------------------------------------------
pimage(m)
pimage(m, order = "AOE")

## -----------------------------------------------------------------------------
pimage(m, order = "AOE", col = rev(bluered()), diag = FALSE, upper_tri = FALSE)
pimage(m, order = "AOE", col = colorRampPalette(c("red", "white", "darkgreen"))(100))

## -----------------------------------------------------------------------------
library("ggplot2")

red_blue <- scale_fill_gradient2(
    low = scales::muted("red"),
    mid = "white",
    high = scales::muted("blue"),
    na.value = "white",
    midpoint = 0)

ggpimage(m, order = "AOE", diag = FALSE, upper_tri = FALSE) + red_blue
  
ggpimage(m, order = "AOE") + scale_fill_gradient2(low = "red", high = "darkgreen")

## -----------------------------------------------------------------------------
d <- as.dist(sqrt(1 - m))

o <- seriate(d, "MDS")
pimage(m , order = c(o, o), main = "MDS", col = rev(bluered()))

o <- seriate(d, "ARSA")
pimage(m , order = c(o, o), main = "ARSA", col = rev(bluered()))

o <- seriate(d, "OLO")
pimage(m , order = c(o, o), main = "OLO", col = rev(bluered()))

o <- seriate(d, "R2E")
pimage(m , order = c(o, o), main = "R2E", col = rev(bluered()))

## -----------------------------------------------------------------------------
if (!require("corrgram")) install.packages("corrgram")
library("corrgram")

corrgram(m, order = "OLO")
corrgram(m, order = "OLO", lower.panel=panel.shade, upper.panel=panel.pie)

## -----------------------------------------------------------------------------
if (!require("corrr")) install.packages("corrr")
library("corrr")

x <- datasets::mtcars |>
       correlate() |>   
       focus(-cyl, -vs, mirror = TRUE) |>  # remove 'cyl' and 'vs'
       rearrange(method = "R2E") |>  
       shave()

rplot(x)

## -----------------------------------------------------------------------------
if (!require("corrplot")) install.packages("corrplot")
library("corrplot")

d <- as.dist(sqrt(1 - m))
o <- seriate(d, "R2E")
m_R2E <- permute(m, c(o,o))

corrplot(m_R2E , order = "original")

