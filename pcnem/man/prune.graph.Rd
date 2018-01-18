\name{prune.graph}
\alias{prune.graph}

\title{Prunes spurious edges in a phenotypic hierarchy}
\description{A heuristic to prune spurious edges in a pehnotypic hierarchy}
\usage{
prune.graph(g,cutIN=NULL,cutOUT=NULL,quant=.95,verbose=TRUE)
}

\arguments{
    \item{g}{an adjacency matrix or a 'graphNEL' object}
    \item{cutIN}{minimum number of missing in-edges required to cut all in-edges. Default}
    \item{cutOUT}{minimum number of missing out-edges required to cut all out-edges}
    \item{quant}{if 'cutIN' or 'cutOUT' are not assigned, a quantile 'quant' of the distribution of missing in- or out-edges for all nodes is used}
    \item{verbose}{Default: TRUE}
}
\details{
\code{prune.graph} provides a heuristic approach to prune surious edges.
\code{prune.graph} compares the input graph to its transitive closure, and counts for each node how many incoming and outgoing edges are missing.
If the number is bigger than a user-defined cutoff, all incoming (outgoing) edges are removed.  
}
\value{
\item{graph}{the pruned phenotypic hierarchy (a 'graphNEL' object)}
\item{removed}{number of removed edges}
\item{missing.in}{number of missing in-edges for each node}
\item{missing.out}{number of missing out-edges for each node}
}

\author{Florian Markowetz}


\examples{
    # a transitively closed core with two spurious edges
    g <- matrix(0,5,5)
    g[1,2] <- 1
    g[2,c(3,4)] <- 1
    g[3,4] <- 1
    g[4,5] <- 1
    dimnames(g) <- list(LETTERS[1:5],LETTERS[1:5])
    g <- as(g,"graphNEL")
    
    # prune graph
    gP <- prune.graph(g)
    
    # plot
    if(require(Rgraphviz)){
      par(mfrow=c(1,2))    
      plot(g,main="two spurious edges")
      plot(gP$graph,main="pruned")
    }
}
\keyword{graphs}
