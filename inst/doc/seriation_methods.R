## -----------------------------------------------------------------------------
library("seriation")
packageVersion('seriation')

## ----results="hide", warning=FALSE, message=FALSE-----------------------------
register_DendSer()
register_optics()
register_smacof()
register_GA()
register_tsne()
register_umap()
register_vegan()

## ----echo = FALSE-------------------------------------------------------------
print_seriation_method_md <- function(x, add_link = TRUE, heading = 4L, ...) {
    cat("\n\n", paste0(strrep("#", heading), collapse = ""), " ", x$name, "\n", sep = "")
  
    cat(x$description, "\n\n")
    
#    cat("* kind:", x$kind, "\n")
    
    optimizes <- x$optimizes
    
    .make_anchor <- function(x) {
      # remove leading numbers
      x <- gsub("^\\d+", "", x)
      # lower case
      tolower(x)
    }
    
    opt_info <- attr(optimizes, "description")
    if(is.na(optimizes)) {
      optimizes <- "N/A"
    } else {
      
      if (add_link) {
        optimizes <- paste0("[", optimizes, "]", "(seriation_criteria.html#", .make_anchor(optimizes) ,")")
      }
   
    }
      
    if(!is.null(opt_info)) {
        optimizes <- paste0(optimizes, " (", opt_info, ")")
    }
     
    cat("* optimizes:", optimizes, "\n")
    cat("* randomized:", x$randomized, "\n")
    
    if (!is.na(x$registered_by)) {
      cat("* registered by: ", x$registered_by, "\n")
    }

  cat("* control parameters: ")
  .print_control(x$control, trim_values = 100)

  cat("\n---\n\n") 
  
  invisible(x)
}

print_method_group <- function(title, desc, methods, l) {
  cat("\n\n###", title, "\n")
  cat(desc, "\n\n")
  available <- vapply(l, function(x) x$name, character(1))
  methods <- intersect(methods, available)
  for (name in methods)
    print_seriation_method_md(l[[which(available == name)[1L]]])
}


.print_control <- function(control,
                           help = TRUE,
                           trim_values = 30L) {
  if (length(control) < 1L) {
    writeLines("no parameters")
  } else{
    contr <- lapply(
      control,
      FUN = function(x)
        strtrim(paste(deparse(x), collapse = ""), trim_values)
    )

    contr <- as.data.frame(t(as.data.frame(contr)))
    colnames(contr) <- c("default")

    contr <- cbind(contr, help = "N/A")
    if (!is.null(attr(control, "help")))
      for (i in seq(nrow(contr))) {
        hlp <- attr(control, "help")[[rownames(contr)[i]]]
        if (!is.null(hlp))
        contr[["help"]][i] <- hlp
      }
    print(knitr::kable(contr))
  }

}

## ----eval=FALSE---------------------------------------------------------------
# seriate(x, method = "Spectral", control = NULL, rep = 1L)

## ----results='asis', echo = FALSE---------------------------------------------
l <- list_seriation_methods(kind = "dist", names_only = FALSE)
available <- vapply(l, function(x) x$name, character(1))

print_method_group("Dendrogram leaf order",
                   "These methods create a dendrogram using hierarchical clustering and then derive the seriation order from the leaf order in the dendrogram. Leaf reordering may be applied.",
  c("DendSer", grep("^DendSer_", available, value = TRUE),
    "GW", grep("^GW_", available, value = TRUE),
    "HC", grep("^HC_", available, value = TRUE),
    "OLO", grep("^OLO_", available, value = TRUE)), l)

print_method_group("Dimensionality reduction",
                   "Find a seriation order by reducing the dimensionality to 1 dimension. This is typically done by minimizing a stress measure or the reconstruction error.",
  c("MDS", "MDS_angle", "isoMDS", "isomap", "monoMDS", "metaMDS",
    "Sammon_mapping", "MDS_smacof", "tsne", "umap"), l)

print_method_group("Optimization",
                   "These methods try to optimize a seriation criterion directly, typically using a heuristic approach.",
  c("ARSA", "Enumerate", "BBURCG", "BBWRCG", "GA", "GSA", "SGLS",
    "QAP_LS", "QAP_2SUM", "QAP_BAR", "QAP_Inertia", "Spectral",
    "Spectral_norm", "TSP"), l)

print_method_group("Other methods",
                   "",
  c("Identity", "optics", "R2E", "Random", "Reverse", "SPIN_NH",
    "SPIN_STS", "VAT"), l)

## ----results='asis', echo = FALSE---------------------------------------------
l <- list_seriation_methods(kind = "matrix", names_only = FALSE)

print_method_group("Seriating rows and columns simultaneously",
                   "Row and column order influence each other.",
  c("BEA", "BEA_TSP", "BK_unconstrained", "CA"), l)

print_method_group("Seriating rows and columns separately using dissimilarities",
                   "",
  "Heatmap", l)

print_method_group("Seriate rows using the data matrix",
                   "These methods need access to the data matrix instead of dissimilarities to reorder objects (rows). The same approach can be applied to columns.",
  c("LLE", "PCA", "PCA_angle", "AOE", "Mean", "Identity", "Reverse",
    "Random"), l)

