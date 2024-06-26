% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AAA_registry_criterion.R
\docType{data}
\name{registry_for_criterion_methods}
\alias{registry_for_criterion_methods}
\alias{registry_criterion}
\alias{list_criterion_methods}
\alias{get_criterion_method}
\alias{set_criterion_method}
\alias{print.criterion_method}
\title{Registry for Criterion Methods}
\format{
An object of class \code{criterion_registry} (inherits from \code{registry}) of length 21.
}
\usage{
registry_criterion

list_criterion_methods(kind, names_only = TRUE)

get_criterion_method(kind, name)

set_criterion_method(
  kind,
  name,
  fun,
  description = NULL,
  merit = NA,
  control = list(),
  verbose = FALSE,
  ...
)

\method{print}{criterion_method}(x, ...)
}
\arguments{
\item{kind}{the data type the method works on. For example, \code{"dist"},
\code{"matrix"} or \code{"array"}.}

\item{names_only}{logical; return only the method name. \code{FALSE} returns
also the method descriptions.}

\item{name}{the name for the method used to refer to the method in the
function \code{\link[=criterion]{criterion()}}.}

\item{fun}{a function containing the method's code.}

\item{description}{a description of the method. For example, a long name.}

\item{merit}{logical; indicating if the criterion measure is a merit
(\code{TRUE}) or a loss (\code{FALSE}) measure.}

\item{control}{a list with control arguments and default values.}

\item{verbose}{logical; print a message when a new method is registered.}

\item{...}{further information that is stored for the method in the
registry.}

\item{x}{an object of class "criterion_method" to be printed.}
}
\value{
\itemize{
\item \code{list_criterion_method()} results is a vector of character strings with the
names of the methods used for \code{criterion()}.
\item \code{get_criterion_method()} returns a given method in form of an object of class
\code{"criterion_method"}.
}
}
\description{
A registry to manage methods used by \code{\link[=criterion]{criterion()}} to calculate a criterion value given data and a
permutation.
}
\details{
All methods below are convenience methods for the registry named
\code{registry_criterion}.

\code{list_criterion_method()} lists all available methods for a given data
type (\code{kind}). The result is a vector of character strings with the
short names of the methods. If \code{kind} is missing, then a list of
methods is returned.

\code{get_criterion_method()} returns information (including the
implementing function) about a given method in form of an object of class
\code{"criterion_method"}.

With \code{set_criterion_method()} new criterion methods can be added by the
user. The implementing function (\code{fun}) needs to have the formal
arguments \verb{x, order, ...}, where \code{x} is the data object, order is
an object of class \link{ser_permutation_vector} and \code{...} can contain
additional information for the method passed on from \code{\link[=criterion]{criterion()}}. The
implementation has to return the criterion value as a scalar.
}
\examples{
## the registry
registry_criterion

# List all criterion calculation methods by type
list_criterion_methods()

# List methods for matrix
list_criterion_methods("matrix")

# get more description
list_criterion_methods("matrix", names_only = FALSE)

# get a specific method
get_criterion_method(kind = "dist", name = "AR_d")

# Define a new method (sum of the diagonal elements)

## 1. implement a function to calculate the measure
criterion_method_matrix_foo <- function(x, order, ...) {
if(!is.null(order)) x <- permute(x,order)
    sum(diag(x))
}

## 2. Register new method
set_criterion_method("matrix", "DiagSum", criterion_method_matrix_foo,
    description = "Calculated the sum of all diagonal entries", merit = FALSE)

list_criterion_methods("matrix")
get_criterion_method("matrix", "DiagSum")

## 3. use all criterion methods (including the new one)
criterion(matrix(1:9, ncol = 3))
}
\seealso{
This registry uses \link[registry:registry]{registry::registry}.

Other criterion: 
\code{\link{criterion}()}
}
\author{
Michael Hahsler
}
\concept{criterion}
\keyword{misc}
