## -----------------------------------------------------------------------------
set.seed(1234)
library("seriation")

data("iris")
x <- as.matrix(iris[sample(nrow(iris), 30), -5])
d <- dist(x)

## -----------------------------------------------------------------------------
register_DendSer()
register_optics()
register_smacof()

## -----------------------------------------------------------------------------
methods <- sort(list_seriation_methods("dist"))
methods <- setdiff(methods, c("BBURCG", "BBWRCG", "Enumerate", "GSA", "SGD", "SGLS"))
methods 

## -----------------------------------------------------------------------------
orders <- list()
criterion <- list()
for (m in methods) {
  cat(m)
  tm <- system.time(orders[[m]] <- seriate(d, method = m))
  criterion[[m]] <- data.frame(time = tm[1]+tm[2], rbind(criterion(d, orders[[m]]))) 
  cat(" took", tm[1]+tm[2], "sec.\n")
}

criterion <- do.call(rbind, criterion)

## -----------------------------------------------------------------------------
orders <- ser_align(orders)
best_to_worse <- order(criterion[["Gradient_weighted"]], decreasing = TRUE)

orders <- orders[best_to_worse]
criterion <- criterion[best_to_worse, ]

## -----------------------------------------------------------------------------
dst <- ser_dist(orders) 
hc <- permute(hclust(dst), order = "OLO", dist = dst)
plot(hc)

## -----------------------------------------------------------------------------
library(DT)
datatable(round(criterion, 2), extensions = "FixedColumns",
    options = list(paging = TRUE, searching = TRUE, info = FALSE,
      sort = TRUE, scrollX = TRUE, fixedColumns = list(leftColumns = 1))) %>%
    formatRound(columns = colnames(criterion) , mark = "", digits=1)

## -----------------------------------------------------------------------------
for (n in names(orders))
  pimage(d, orders[[n]], main = n , key = FALSE)

## -----------------------------------------------------------------------------
methods <- sort(list_seriation_methods("matrix"))

# AOE if for correlation matrices only
methods <- setdiff(methods, c("AOE"))
methods 

## ----fig.height= 5------------------------------------------------------------
orders <- list()
criterion <- list()
for (m in methods) {
  cat(m)
  tm <- system.time(orders[[m]] <- seriate(x, method = m))
  criterion[[m]] <- data.frame(time = tm[1]+tm[2], rbind(criterion(x, orders[[m]])))
  cat(" took", tm[1]+tm[2], "sec.\n")
}

criterion <- do.call(rbind, criterion)

## -----------------------------------------------------------------------------
datatable(round(criterion, 2), extensions = "FixedColumns",
    options = list(paging = TRUE, searching = TRUE, info = FALSE,
      sort = TRUE, scrollX = TRUE, fixedColumns = list(leftColumns = 1))) %>%
    formatRound(columns = colnames(criterion) , mark = "", digits = 1)

## -----------------------------------------------------------------------------
best_to_worse <- order(criterion[["Moore_stress"]], decreasing = FALSE)

orders <- orders[best_to_worse]
criterion <- criterion[best_to_worse, ]

## -----------------------------------------------------------------------------
for (n in names(orders))
  pimage(x, orders[[n]], main = n , key = FALSE)

