\name{SCCgraph}
\alias{SCCgraph}

\title{Combines Strongly Connected Components into single nodes}
\description{
  \code{SCCgraph} is used to identify all nodes which are not distinguishable given the data. 
}
\usage{
SCCgraph(x,name=TRUE,nlength=20)
}

\arguments{
  \item{x}{graphNEL object or an adjacency matrix}
  \item{name}{Concatenate all names of summarized nodes, if TRUE, or number nodes, if FALSE. Default: TRUE}
  \item{nlength}{maximum length of names}
}
\details{
  A graph inferred by either \code{nem} or \code{nemModelSelection} may have cycles if some phenotypic profiles are not distinguishable.
  The function \code{SCCgraph} identifies cycles in the graph (the strongly conneced components) and summarizes them in a single node.
  The resulting graph is then acyclic.  
}

\value{
 \item{graph}{a graphNEL object with connected components of the input graph summarized into single nodes}
 \item{scc}{a list mapping SCCs to nodes}
 \item{which.scc}{a vector mapping nodes to SCCs}
}

\author{Florian Markowetz, Holger Froehlich}

\seealso{\code{\link{nem}}, \code{\link{transitive.reduction}}}
\examples{
   data("BoutrosRNAi2002")
   D   <- BoutrosRNAiDiscrete[,9:16]
   res <- nem(D,control=set.default.parameters(unique(colnames(D)), para=c(.13,.05)))
   # 
   sccg <- SCCgraph(res$graph,name=TRUE)
   #
   par(mfrow=c(1,2))
   if(require(Rgraphviz)){
    plot.nem(res, main="inferred from data")      
    plot(sccg$graph, main="condensed (rel,key)")
    }
}
\keyword{graphs}
