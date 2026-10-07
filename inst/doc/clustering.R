## ----message=FALSE, warning=FALSE---------------------------------------------
library("seriation")
library("ggplot2")
library("dplyr")

set.seed(1234)

## ----ruspini------------------------------------------------------------------
library(cluster)
data(ruspini)
ruspini <- ruspini |> sample_frac()
head(ruspini)

plot(ruspini)

## ----ruspini4-----------------------------------------------------------------
cl_ruspini <- kmeans(ruspini, centers = 4, nstart = 5)

d_ruspini <- ruspini |> dist()
ggdissplot(d_ruspini, cl_ruspini$cluster) + ggtitle("Dissimilarity Plot")

clusplot(ruspini, cl_ruspini$cluster, labels = 4)

## ----ruspini3-----------------------------------------------------------------
cl_ruspini3 <- kmeans(ruspini, center=3, nstart = 5)

ggdissplot(d_ruspini, cl_ruspini3$cluster) + ggtitle("Dissimilarity Plot")

clusplot(ruspini, cl_ruspini3$cluster, labels = 4)

## ----ruspini7-----------------------------------------------------------------
cl_ruspini7 <- kmeans(ruspini, centers=7)

ggdissplot(d_ruspini, cl_ruspini7$cluster) + ggtitle("Dissimilarity Plot")

clusplot(ruspini, cl_ruspini7$cluster)

## ----ruspini0-----------------------------------------------------------------
ggdissplot(d_ruspini) + ggtitle("Dissimilarity Plot")

## ----message=FALSE------------------------------------------------------------
library(cluster)

data(Votes, package = "cba")
x <- cba::as.dummy(Votes[-17])
d_votes <- dist(x, method = "binary")

## ----votes2-------------------------------------------------------------------
labels_votes2 <- pam(d_votes, k=2, cluster.only = TRUE)

ggdissplot(d_votes, labels_votes2) + ggtitle("Dissimilarity Plot")

clusplot(d_votes, diss = TRUE, labels_votes2, labels = 4)

## ----votes12------------------------------------------------------------------
labels_votes12 <- pam(d_votes, k=12, cluster.only = TRUE)

ggdissplot(d_votes, labels_votes12) + ggtitle("Dissimilarity Plot")

clusplot(d_votes, diss = TRUE, labels_votes12, labels = 4)

