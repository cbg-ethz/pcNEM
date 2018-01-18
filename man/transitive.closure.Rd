\name{transitive.closure}
\alias{transitive.closure}

\title{Computes the transitive closure of a directed graph}
\description{
  Computes the transitive closure of a graph. Introduces a direct edge whenever 
  there is a path between two nodes in a digraph.
}
\usage{
transitive.closure(g, mat=FALSE, loops=TRUE)
}

\arguments{
  \item{g}{graphNEL object or adjacency matrix.}
  \item{mat}{convert result to adjacency matrix.}
  \item{loops}{Add loops from each node to itself?}
}
\details{
  This function calculates the transitive closure of a given
  graph. We use the matrix exponential to find the transitive closure.
}
\value{
returns a graphNEL object or adjacency matrix
}

\author{Florian Markowetz}

\seealso{\code{\link{transitive.reduction}}}
\examples{
   V <- LETTERS[1:3]
   edL <- list(A=list(edges="B"),B=list(edges="C"),C=list(edges=NULL))
   g <- new("graphNEL",nodes=V,edgeL=edL,edgemode="directed")
   gc <- transitive.closure(g,loops=FALSE)
    
  if(require(Rgraphviz)){
    par(mfrow=c(1,2))
    plot(g,main="NOT transitively closed")
    plot(gc,main="transitively closed")
  }
}
\keyword{graphs}
