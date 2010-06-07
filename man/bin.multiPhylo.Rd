\name{bin.multiPhylo}
\Rdversion{1.1}
\alias{bin.multiPhylo}
\title{
    Bin Trees
}

\description{
    Bins trees according to branching topology. 
}

\usage{
    bin.multiPhylo(treelist)
}

\arguments{
  \item{treelist}{
      A list of trees that can be passed to dist.phylo (see the help for
      dist.phylo for acceptable formats)
}
}

\details{
   Bins trees according to branching topology. Two trees are considered to have
   the same topology if the same set of partitions of tips are produced by
   the edges, which corresponds to the same branching up to rearrangement
   of tips.
}
\value{
    Returns a numeric vector of bin ids. Bin ids are assigned in order of the
    first tree in that bin, that is, the first k unique trees in the list
    passed will be assigned bins 1..k in order of appearance. 
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
breps = 500

trees <- boot.phylo(otree, woodmouse, B=breps, function(x)
        root(fastme.ols(dist.dna(x)), "No305", resolve.root=TRUE),trees=TRUE)

combined.trees <- c(list(otree), trees$trees)

bin.multiPhylo(combined.trees)

}
