\name{transitive.projections}
\alias{transitive.projections}
\alias{CheckEdge}
\alias{remTwoEdges}
\alias{transSubGr}
\alias{EdgeEk}
\alias{distdecrease}
\alias{distsame}
\alias{distincrease}
\alias{distincrease1}
\alias{OneNeighborhood}
\alias{TwoNeighborhood}
\alias{ThreeNeighborhood}
\alias{FourNeighborhood}
\alias{VecToMat}
\alias{is.transitive}

\title{Computes the transitive approximation of a directed graph}
\description{
  Computes the transitive approximation of a graph. The transitive approximation of a graph is a graph that is "almost" transitively closed and has minimal distance to the input graph.
}
\usage{
transitive.projections(adjmat)
}

\arguments{
  \item{adjmat}{graphNEL object or adjacency matrix.}
}
\value{
returns adjacency matrices and having minimal graph distance to the input graph matrix
}

\author{Juby Jacob}

\seealso{\code{\link{transitive.projections}}}
\keyword{graphs}
