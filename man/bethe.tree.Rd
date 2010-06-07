\name{bethe.tree}
\Rdversion{1.1}
\alias{bethe.tree}
\title{
    Bethe Tree
}

\description{
    Generates a Bethe tree with given tips, inner edge lengths, and outgroup.
}

\usage{
    bethe.tree(tips, level.lengths = NULL, outgroup="O", outgroup.dist=1)
}

\arguments{
  \item{tips}{
      A list of tip names as a character vector. Should be a power of 2. All
      tip names must be distinct.
}
 \item{level.lengths}{
     Edge lengths for each level, counted from the bottom up. NULL means a
     default of 1. If the vector isn't long enough, the last value will be
     repeated as necessary.
}
 \item{outgroup}{
     The tip label for the outgroup.
}
 \item{outgroup.dist}{
     The distance of the outgroup from the root. 
}
}

\details{
    Generates a Bethe tree with specified internal edge lengths.
}
\value{
    A class of type \code{phylo} representing the tree.
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

plot(bethe.tree( as.character(1:16), 1:4, "17", 14)) 

}
