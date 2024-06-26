% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_order.R
\name{get_order}
\alias{get_order}
\alias{get_order.ser_permutation_vector}
\alias{get_order.ser_permutation}
\alias{get_order.hclust}
\alias{get_order.dendrogram}
\alias{get_order.integer}
\alias{get_order.numeric}
\alias{get_rank}
\alias{get_permutation_matrix}
\title{Extracting Order Information from a Permutation Object}
\usage{
get_order(x, ...)

\method{get_order}{ser_permutation_vector}(x, ...)

\method{get_order}{ser_permutation}(x, dim = 1, ...)

\method{get_order}{hclust}(x, ...)

\method{get_order}{dendrogram}(x, ...)

\method{get_order}{integer}(x, ...)

\method{get_order}{numeric}(x, ...)

get_rank(x, ...)

get_permutation_matrix(x, ...)
}
\arguments{
\item{x}{an object of class \link{ser_permutation} or
\link{ser_permutation_vector}.}

\item{...}{further arguments are ignored for \code{get_order()}.  For
\code{get_rank()} and for \code{get_permutation_matrix()} the additional
arguments are passed on to \code{get_order()} (e.g., as \code{dim}).}

\item{dim}{order information for which dimension should be returned?}
}
\value{
Returns an integer permutation vector/a permutation matrix.
}
\description{
Method to get the order information from an object of class
\link{ser_permutation} or \link{ser_permutation_vector}. Order information
can be extracted as a permutation vector, a vector containing each
object's rank or a permutation matrix.
}
\details{
\code{get_order()} returns the permutation as an integer vector which arranges the
objects in the seriation order. That is, a vector with the index of the first,
second, \eqn{..., n}-th object in the order defined by the permutation.
These permutation vectors can directly be
used to reorder objects using subsetting with \code{"["}.  \emph{Note:} In
\pkg{seriation} we usually use these order-based permutation vectors.
\strong{Note on names:} While R's \code{\link[=order]{order()}} returns an unnamed vector,
\code{get_order()} returns names (if available). The names are the object label
corresponding to the index at that position.
Therefore, the names in the order are in the order after
the permutation.

\code{get_rank()} returns the seriation as an integer vector containing the
rank/position for each objects after the permutation is applied.
That is, a vector with the position of the first, second,
\eqn{..., n}-th object after permutation.  Note: Use
\code{order()} to convert ranks back to an order.

\code{get_permutation_matrix()} returns a \eqn{n \times n}{n x n} permutation
matrix.
}
\examples{
## create a random ser_permutation_vector
## Note that ser_permutation_vector is a single permutation vector
x <- structure(1:10, names = paste0("X", 1:10))
o <- sample(x)
o

p <- ser_permutation_vector(o)
p

get_order(p)
get_rank(p)
get_permutation_matrix(p)

## reorder objects using subsetting, the provided permute function or by
## multiplying the with the permutation matrix. We use here
x[get_order(p)]
permute(x, p)
drop(get_permutation_matrix(p) \%*\%  x)

## ser_permutation contains one permutation vector for each dimension
p2 <- ser_permutation(p, sample(5))
p2

get_order(p2, dim = 2)
get_rank(p2, dim = 2)
get_permutation_matrix(p2, dim = 2)
}
\seealso{
Other permutation: 
\code{\link{permutation_vector2matrix}()},
\code{\link{permute}()},
\code{\link{ser_dist}()},
\code{\link{ser_permutation}()},
\code{\link{ser_permutation_vector}()}
}
\author{
Michael Hahsler
}
\concept{permutation}
\keyword{manip}
