\name{distory-package}
\Rdversion{1.1}
\alias{distory-package}
\alias{distory}
\docType{package}
\title{
    Distance Between Phylogenetic Histories
}
\description{
    The 'distory' package provides functions for computing geodesic distances
    between phylogenetic trees, as well as functions which use this distance.
    Methods for computing Gromov delta-hyperbolicity, Markov Chain Monte Carlo
    routines in tree space, and per-position leverage for DNA sequences are
    included.
}
\details{
    A description of the algorithm used for the distance computation can be found
    in help(dist.multiPhylo).
}
\author{
    John Chakerian <chakj@stanford.edu> and Susan Holmes <susan@stat.stanford.edu>
}
\references{
    Chakerian, J. and Holmes, S. P. Computational Tools for Evaluating
    Phylogenetic and Heirarchical Clustering Trees.	arXiv:1006.1015v1. 

    Billera, L. J., Holmes, S. P., and Vogtmann, K. (2001) Geometry of the space
    of phylogenetic trees. _Adv. Appl. Math_, *27*, 733-767.

    Megan Owen, J. Scott Provan, "A Fast Algorithm for Computing Geodesic
    Distances in Tree Space," IEEE/ACM Transactions on Computational Biology
    and Bioinformatics, 14 Jan. 2010.
}
\keyword{ package }
\examples{

data(woodmouse)
otree <- root(nj(dist.dna(woodmouse)), "No305", resolve.root=TRUE)
breps = 250

trees <- boot.phylo(otree, woodmouse, B=breps, function(x)
        root(nj(dist.dna(x)), "No305", resolve.root=TRUE),trees=TRUE)

combined.trees <- c(list(otree), trees$trees)
tree.dists <- dist.multiPhylo(combined.trees)

mdres <- cmdscale(tree.dists, k=breps, add=TRUE)
plot(mdres$points[,1], mdres$points[,2], col = c("red", rep("black", breps)))
text(mdres$points[,1], mdres$points[,2], labels=1:(breps+1), cex=0.7, adj=c(0,2))

}
