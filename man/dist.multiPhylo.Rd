\name{dist.multiPhylo}
\Rdversion{1.1}
\alias{dist.multiPhylo}
\title{
    Geodesic Distance Between Phylogenetic Trees
}

\description{
    Computes the geodesic distance of a list of phylogenetic trees using a 
    polynomial algorithm. 
}

\usage{
    dist.multiPhylo(x, method="geodesic", force.multi2di = FALSE, outgroup = NULL, 
    convert.multifurcating = FALSE, use.random.resolution = FALSE, scale =
    NULL, verbose = FALSE)
}

\arguments{
  \item{x}{
      A list of ape trees (class 'phylo'). The list does not have to be of 
      class 'multiPhylo'. The function will also accept a list of strings of 
      trees in Newick format, or a single string with trees in Newick format
      separated by semicolons. All the trees must have the same tip labels.
}
  \item{method}{
      Determines which distance method is used. Options are 'geodesic' for the
      tree space geodesic distance, or 'edgeset' for the number of edges
      (defined by splits of tips) that are different.
}
  \item{force.multi2di}{
      Force conversion of every tree to strict bifurcating through the ape 
      function 'multi2di', using the use.random.resolution as its parameter.
      This option should not be used in conjunction with specification of an
      outgroup.
}
  \item{outgroup}{
      Specifies an outgroup to root each tree with respect to. This calls the
      ape function 'root' on every tree in the list.
}
  \item{convert.multifurcating}{
      Setting this option will check every tree for multifurcations using the 
      ape function 'is.binary.tree' - if it returns FALSE, the ape function
      'multi2di' will be called on it. Note that this does not ensure a tree
      is strictly binary, since ape considers an unrooted tree binary even if
      the root node is trifurcating. This option can be used in conjunction
      with specification of an outgroup.
}

  \item{use.random.resolution}{
      Specifies the parameter to 'multi2di' if needed. 
}
  \item{scale}{
      Specifies a scale to make all trees unformly scaled (that is, the sum of
      all edges will be uniform)scale to make all trees unformly scaled (that 
      is, the sum of all edge lengths will be uniform). The parameter can 
      either be a tree of class \code{phylo} or a numeric value for the sum of 
      all edge lengths.
}
  \item{verbose}{
      Turns on incremental status updates and more warnings. Helpful for large 
      computations.
}

}
\details{
    This function computes the geodesic distance according to Billera et. al.
    using an algorithm based off of the polynomial time algorithm of Owen
    and Provan. Since it corresponds to a formal definition of tree-space as
    a space of strictly binary trees, no mulifurcations are allowed, including
    on the root node. In addition, negative and 0-lengthed edges are clamped to
    a very small value (DBL_MIN) for technical reasons. 

    The Newick parser supports only a subset of the Newick format. In
    particular, it does not at the moment allow for internal node labels, only
    weights. Weights will be automatically set to 1 if not specified. It may be
    necessary to clean data in ape to make the trees conform to this.
}
\value{
    Returns a distance matrix of class 'dist' representing the pairwise
    geodesic distances between all input trees. Keep in mind this distance
    matrix is not Euclidean. N/A values are provided in the case of an error in
    determining the distance.
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
\author{
    John Chakerian <chakj@stanford.edu>
}

\seealso{
    \code{\link[ape]{dist.dna}}
    \code{\link[ape]{boot.phylo}}
    \code{\link{cmdscale}} 
}

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
text(mdres$points[,1], mdres$points[,2], ,labels=1:(breps+1), cex=0.7, adj=c(0,2))

}
