\name{closest.transitive.greedy}
\alias{closest.transitive.greedy}

\title{Find transitively closed graph most similar to the given one}
\description{
	First, from the original graph \eqn{\Phi} spurious edges are pruned via
	\code{prune.graph}. Then the new graph \eqn{\Phi'} is transitively
	closed. Afterwards, the algorithms successively introduces new edges
	minimizing the distance to the original graph (defined as
	\eqn{\sum_{ij} |\Phi_{ij} - \Phi'_{ij}|}) most. After each edge
	addition the graph is transitively closed again. 
}
\usage{
closest.transitive.greedy(Phi, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{  
\item{Phi}{adjacency matrix}    
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}  
}
\value{
  adjacency matrix
}
\author{Holger Froehlich}

\seealso{\code{\link{prune.graph}}, \code{\link{transitive.closure}}, \code{\link{transitive.reduction}}}

\keyword{models}
