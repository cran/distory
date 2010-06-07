\name{mcmc.target.seq}
\Rdversion{1.1}
\alias{mcmc.target.seq}
\alias{boot.samples.idxs}
\alias{lookup.samples}
\alias{convert.table.to.idx}

\title{
    Find MCMC Target Sequence
}

\description{
    \code{mcmc.target.seq} uses MCMC to find a configuration of DNA positions
    to get as close as possible to a given tree.

    \code{boot.samples.idxs} bootstraps over indices into a DNA matrix.

    \code{lookup.samples} goes from an index representation of a configuration
    of DNA to the actual DNAbin format.

    \code{convert.table.to.idx} converts a table of counts for positions 1..n
    into a list of indices corresponding to positions (i.e. goes from the
    tabled form to a vector whose tabling matches the input).
}

\usage{
    mcmc.target.seq(data, x, F, n) 

    boot.samples.idxs(data, B = 100, block = 1)

    lookup.samples(data, idxs)

    convert.table.to.idx(T)
}

\arguments{
  \item{data}{
      A DNA matrix in DNAbin format.
}
  \item{x}{
      A tree of class 'phylo' to estimate.
}
  \item{F}{
      A tree estimation function, accepting a DNA matrix in DNAbin format and
      returning a tree of class 'phylo.'
}
  \item{n}{
      The number of MCMC iterations to perform
}
 \item{B}{
      The number of bootstrap replicates.
}
 \item{block}{
      The block size to use during bootstrapping.
}
 \item{idxs}{
     A list of numeric vectors of indices to use for lookup.
}
 \item{T}{
     A table or table-like vector to convert.
}
}

\details{
    \code{mcmc.target.seq} performs an MCMC with simulated annealing to locate
    a configuration of DNA positions from the original matrix that gets as
    close as possible to a target tree. Propositions for the MCMC replacing
    one character with another uniformly at random. 

    The remaining functions are intended to be used as support functions.
}
\value{
    \code{mcmc.target.seq} returns a list of 4 elements: a numeric vector of
    counts of each position in the original matrix, the best estimated
    tree, the best distance from the estimated tree to the target tree, and
    a numeric vector of the distances for every iteration of the
    simulation. 

    \code{boot.samples.idxs} returns a numeric vector representing the
    bootstrapped idices.

    \code{lookup.samples} returns a list of objects of class DNAbin
    corresponding to the DNA sequences generated from indices into the original
    DNA matrix.

    \code{convert.table.to.idx} returns a numeric vector of indices based on
    the table counts.
}
\author{
    John Chakerian <chakj@stanford.edu>
}
\references{
    Chakerian, J. and Holmes, S. P. Computational Tools for Evaluating
    Phylogenetic and Heirarchical Clustering Trees.	arXiv:1006.1015v1.
}

\seealso{
    \code{\link[distory]{dist.multiPhylo}}
    \code{\link[distory]{orthant.boundary.tree}}
}

\examples{

data(woodmouse)
otree <- root(fastme.ols(dist.dna(woodmouse)), "No305", resolve.root=TRUE)
breps = 200

trees <- boot.phylo(otree, woodmouse, B=breps, function(x)
        root(fastme.ols(dist.dna(x)), "No305", resolve.root=TRUE),trees=TRUE)

combined.trees <- c(list(otree), trees$trees)

binning <- bin.multiPhylo(combined.trees)

tree.a <- combined.trees[[match(1, binning)]]
i = 2
max.bin = max(binning)
tree.b <- combined.trees[[match(2, binning)]]

while(length(distinct.edges(tree.a,tree.b)) > 1 && i < max.bin)
{
    i = i + 1
    tree.b = combined.trees[[match(i, binning)]]
}

bdy.tree <- orthant.boundary.tree(tree.a, tree.b)

f.est <- function(x) root(nj(dist.dna(x)), "No305", resolve.root=TRUE)

res <- mcmc.target.seq(woodmouse, bdy.tree, f.est, 1000)

par(mfrow=c(2,1))
plot(res$tree)
plot(res$vals)

}
