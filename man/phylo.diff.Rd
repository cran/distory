\name{phylo.diff}
\Rdversion{1.1}
\alias{phylo.diff}
\alias{distinct.edges}
\alias{edge.from.split}
\alias{get.bipartition}
\alias{partition.leaves}

\title{
    Differences Between Phylogenetic Trees
}

\description{
    A family of functions for determining and plotting the differences between
    two trees.

    \code{phylo.diff} plots two trees side by side, highlighting edges unique
    to each tree in red.

    \code{distinct.edges} finds the edges present in the first argument not in
    the second.

    \code{edge.from.split} locates the edge id from a given split.

    \code{get.bipartition} gets the bipartition of tips formed by a single edge.

    \code{partition.leaves} returns the set of all bipartitions from all edges.
}

\usage{
    phylo.diff(x, y, ...)

    distinct.edges(x, y)

    edge.from.split(x, split)

    get.bipartition(x, e)

    partition.leaves(x)

}

\arguments{
  \item{x}{
      The first (or only) tree.
}
  \item{y}{
      The second tree, for the functions that accept two trees.
}
  \item{split}{
      A list of bipartitions, probably from \code{partition.leaves}.
}
  \item{e}{
      An edge for a particular tree, given as an id.
}
 \item{...}{
      Additional arguments to pass to the \code{plot.phylo} function.
}
}

\details{
    \code{phylo.diff} uses the ape tree plotting function. The other functions
    are mostly meant as support functions.
}
\value{
    \code{phylo.diff} returns invisible.

    \code{distinct.edges} returns a numeric vector of edge ids for the first
    tree.

    \code{edge.from.split} returns an edge id for a particular tree
    corresponding to a given bipartition and NA if none such edge exists.

    \code{get.bipartition} returns a character vector of the tips below that
    edge in the given tree.

    \code{partition.leaves} returns a list of partitions (themselves character
    vectors) of the given tree.
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
}

\examples{

data(woodmouse)
otree <- root(fastme.ols(dist.dna(woodmouse)), "No305", resolve.root=TRUE)
breps = 10

trees <- boot.phylo(otree, woodmouse, B=breps, function(x)
        root(fastme.ols(dist.dna(x)), "No305", resolve.root=TRUE),trees=TRUE)

combined.trees <- c(list(otree), trees$trees)

binning <- bin.multiPhylo(combined.trees)

phylo.diff(combined.trees[[match(1, binning)]], combined.trees[[match(2, binning)]])

}
