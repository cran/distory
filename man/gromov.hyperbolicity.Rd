\name{gromov.hyperbolicity}
\Rdversion{1.1}
\alias{gromov.hyperbolicity}
\title{
    Gromov Hyperbolicity Constant
}

\description{
    Computes the Gromov Hyperbolicity Constant of a distance matrix.
}

\usage{
    gromov.hyperbolicity(d, deltas = FALSE, scale = NA)
}

\arguments{
  \item{d}{
      A distance matrix of type \code{dist} or \code{matrix}, or anything that
      can be coerced into \code{dist} by \code{as.dist}. Must have at least 4
      points.
}
  \item{deltas}{
      A logical value specifying whether to return the vector of delta values.
      Default is \code{FALSE}.
}
  \item{scale}{
     Specifies a scaling method for each delta. Default is no scaling (NA or
     "none"). Available methods are "max" which scales deltas by the max of the
     sums computed, and "perimeter" for the largest perimeter of the four points.
}
}

\details{
    This computes a constant that represents the relaxation of a 4-point
    condition for delta-hyperbolicity. See (Gromov 1987) for details.
}

\value{
    The Gromov hyperbolicity constant of the given distance matrix.
}

\author{
    John Chakerian <chakj@stanford.edu>
}

\seealso{
    \code{\link[distory]{dist.multiPhylo}}
}

\references{
    M. Gromov. Hyperbolic groups. In Essays in Group Theory, pages 73-263.
    Springer, New York, 1987.

    Chakerian, J. and Holmes, S. P. Computational Tools for Evaluating
    Phylogenetic and Heirarchical Clustering Trees.	arXiv:1006.1015v1. 
}

\examples{

# scale final delta by max distance
points <- cbind(runif(100), runif(100))
d <- dist(points)
gromov.hyperbolicity(d)/max(d)

# scale each delta by max distance for the 4 points
points <- cbind(runif(100), runif(100))
d <- dist(points)
gromov.hyperbolicity(d, scale="max")

# scale each delta by the max perimeter for the 4 points
points <- cbind(runif(100), runif(100))
d <- dist(points)
gromov.hyperbolicity(d, scale="max")

}
