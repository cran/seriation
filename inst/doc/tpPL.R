## -----------------------------------------------------------------------------
library("seriation")
data("Chameleon")

x <- chameleon_ds7[sample(1:nrow(chameleon_ds7), 500), ]
plot(x)

## -----------------------------------------------------------------------------
D <- dist(x)

hc <- hclust(D, method = "complete")
plot(hc, labels = FALSE)

P <- cophenetic(hc)

## -----------------------------------------------------------------------------
beta <- 1
F <- D + beta * P

o <- seriate(F, method = "TSP")
criterion(F, o, method = "Path_length")

## -----------------------------------------------------------------------------
betas <- c(0, .25, .5, 1, 2, 10)

for (beta in betas) {
  F <- D + beta * P
  pimage(F, order = "TSP", main = paste("beta =", beta))
}

