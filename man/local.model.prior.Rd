\name{local.model.prior}
\alias{local.model.prior}
\title{Computes a prior to be used for edge-wise model inference}
\description{
  The function \code{pairwise.posterior} infers a phenotypic hierarchy edge by edge by 
  choosing between four models (unconnected, subset, superset, undistinguishable). 
  For each edge, \code{local.model.prior} computes a prior distribution over the four models.
  It can be used to ensure sparsity of the graph and high confidence in results.
}
\usage{
local.model.prior(size,n,bias)
}

\arguments{
  \item{size}{expected number of edges in the graph. }
  \item{n}{number of perturbed genes in the dataset, number of nodes in the graph}
  \item{bias}{the factor by which the double-headed edge is preferred over the single-headed edges}
}
\details{
 A graph on \code{n} nodes has \code{N=n*(n-1)/2} possible directed edges (one- or bi-directional). 
 If each edge occurs with probability \eqn{p}, we expect to see \eqn{Np} edges in the graph.
 The function \code{local.model.prior} takes the number of genes (\code{n}) and the 
 expected number of edges (\code{size}) as an input and computes a prior distribution 
 for edge occurrence: no edge with probability \code{size/N}, and the probability for 
 edge existence being split over the three edge models with a bias towards the conservative 
 double-headed model specified by \code{bias}. To ensure sparsity, the \code{size} should 
 be chosen small compared to the number of possible edges.
}
\value{
    a distribution over four states: a vector of four positive real numbers summing to one
}

%\references{}
\author{Florian Markowetz}
%\note{}

\seealso{\code{\link{pairwise.posterior}}, \code{\link{nem}}}
\examples{
# uniform over the 3 edge models
local.model.prior(4,4,1)
# bias towards <->
local.model.prior(4,4,2)
}
\keyword{models}
