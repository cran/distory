\name{position.leverage}
\Rdversion{1.1}
\alias{position.leverage}
\title{Position Leverage}
\description{
  Provides a rough heuristic for determining the degree to which each
  position in the DNA matrix affects the resulting tree.
}
\usage{
position.leverage(data, F, to = NULL, rep = 50, by = 1)
}
\arguments{
  \item{data}{A DNA matrix in \code{DNAbin} format.}
  \item{F}{
    A tree estimation function, accepting a DNA matrix of class \code{DNAbin}
    and returning a tree of class \code{phylo}.}
  \item{to}{
    The tree with which distances are measured in respect to, or \code{NULL}
    to indicate the tree estimated by \code{F} for the starting DNA matrix.}
  \item{rep}{
    The number of times to replicate the position in question.}
  \item{by}{
    The function will perform the calculation on every \code{by}-th position
    (that is, on \code{seq(1, N, by)}).}}
\details{
  This function takes a DNA matrix and, for every \code{by}-th position,
  replicates that position \code{rep} times, randomly removing
  \code{rep} other positions to keep all sequences the same length
  other positions to keep all sequences the same length. For each new DNA
  matrix created in this way, \code{F} is used to estimate the corresponding
  tree, and the distance to tree \code{to} is computed and stored. This
  distance can be thought of as somewhat analogous to the leverage of that
  position.
}
\value{
  Returns a numeric vector of distances from tree \code{to} for each position
  sampled.
}
\author{John Chakerian}
\references{
  Chakerian, J. and Holmes, S. P. Computational Tools for Evaluating
  Phylogenetic and Heirarchical Clustering Trees. arXiv:1006.1015v1.
}
\seealso{\code{\link{dist.multiPhylo}}}
\examples{
data(woodmouse)
f.est <- function(x) root(nj(dist.dna(x)), "No305", resolve.root = TRUE)
position.leverage(woodmouse, f.est, by = 10)
}
\keyword{manip}
