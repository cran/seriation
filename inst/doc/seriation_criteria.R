## -----------------------------------------------------------------------------
library("seriation")
packageVersion('seriation')

## ----echo = FALSE-------------------------------------------------------------
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
    colnames(contr) <- "default"

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

optimized_by <- function(criterion, add_link = TRUE) {
  l <- registry_seriate$get_entries()
  names <- sapply(l, "[[", "name")
  ops <- sapply(l, "[[", "optimizes")
  match <- ops == criterion
  match[is.na(match)] <- FALSE
  res <- unname(names[match])
  if (length(res) == 0) 
    return("N/A")
  
  .make_anchor <- function(x) {
    # remove leading numbers
    x <- gsub("^\\d+", "", x)
    # lower case
    tolower(x)
  }
  
  # add links 
  if (add_link)
    res <- sapply(res, FUN = function(m) paste0("[", m , "]", "(seriation_methods.html#", .make_anchor(m), ")"))
  
  res
}

print_criterion_method_md <- function(x, ...) {
    cat("\n\n###", x$name, "\n")
    cat("\n", x$description, "\n\n")
  
#    cat("* kind:", x$kind, "\n")
    cat("* merit:", x$merit, "\n")
    cat("* optimized by: ", paste(optimized_by(x$name), collapse = ", "), "\n")
    if (!is.na(x$registered_by)) {
      cat("* registered by: ", x$registered_by, "\n")
    }
    
    #    cat("* randomized:", x$randomized, "\n\n")
    

  writeLines("* additional parameters:")
  .print_control(x$control, trim_values = 100)

  cat("\n---\n\n")
  
  invisible(x)
}


## ----results="hide", warning=FALSE, message=FALSE-----------------------------
register_DendSer()
register_smacof()

## ----results='asis', echo = FALSE---------------------------------------------
l <- list_criterion_methods(kind = "dist", names_only = FALSE)
l <- l[order(sapply(l, "[[", "name"))]
for(i in seq_along(l))
  print_criterion_method_md(l[[i]])

## ----results='asis', echo = FALSE---------------------------------------------
l <- list_criterion_methods(kind = "matrix", names_only = FALSE)
l <- l[order(sapply(l, "[[", "name"))]
for(i in seq_along(l))
  print_criterion_method_md(l[[i]])

