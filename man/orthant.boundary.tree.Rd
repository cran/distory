\name{orthant.boundary.tree}
\Rdversion{1.1}
\alias{orthant.boundary.tree}
\title{
    Orthant Boundary Tree
}

\description{
    Produces a degenerate tree on the boundary between trees that differ by one split.
}

\usage{
    orthant.boundary.tree(x,y)
}

\arguments{
  \item{x}{
      The tree in the first orthant.
}
  \item{y}{
      The tree in the second orthant.
}
}

\details{
    The tree found is the tree on the boundary between the two orthants such
    that it is on the straight line connecting the two trees when one
    orthant is thought of as being the (-,+) quadrant and the second
    orthant as being the (+,+) quadrant, where the (0,y) line is the particular
    boundary in question.
}
\value{
    Returns an object of class 'phylo' representing the boundary tree.
}
\author{
    John Chakerian <chakj@stanford.edu>
}
\references{
    Chakerian, J. and Holmes, S. P. Computational Tools for Evaluating
    Phylogenetic and Heirarchical Clustering Trees.	arXiv:1006.1015v1.
}

\seealso{
    \code{\link[distory]{mcmc.target.seq}}
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

plot(orthant.boundary.tree(tree.a, tree.b))

}
