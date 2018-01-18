\name{transitive.reduction}
\alias{transitive.reduction}

\title{Computes the transitive reduction of a graph}
\description{
  \code{transitive.reduction} removes direct edges, which can be explained by another path in the graph. Regulation directions  inferred via \code{infer.edge.type} are taken into account.
}
\usage{
transitive.reduction(g)
}

\arguments{
  \item{g}{adjacency matrix}
}
\details{
  \code{transitive.reduction} uses a modification of the classical algorithm from the Sedgewick book for computing transitive closures. The so-called "transitive reduction" is neither necessarily unique (only for DAGs) nor minimal in the number of edges (this could be improved). 
}
\value{
returns an adjacency matrix with shortcuts removed
}
\references{
R. Sedgewick, Algorithms, Pearson, 2002.
}

\author{Holger Froehlich}

\seealso{\code{\link{transitive.closure}}, \code{\link{infer.edge.type}}}
\examples{
   V <- LETTERS[1:3]
   edL <- list(A=list(edges=c("B","C")),B=list(edges="C"),C=list(edges=NULL))
   gc <- new("graphNEL",nodes=V,edgeL=edL,edgemode="directed")
   g <- transitive.reduction(gc)
    
   if(require(Rgraphviz)){
    par(mfrow=c(1,2))
    plot(gc,main="shortcut A->C")
    plot(as(g,"graphNEL"),main="shortcut removed")
  }
}
\keyword{graphs}
