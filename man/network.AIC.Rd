\name{network.AIC}
\alias{network.AIC}

\title{AIC/BIC criterion for network graph}
\description{
  Calclate AIC/BIC for a given network graph (should be transitively closed). The number of free parameters equals the number of unknown edges in the network graph.
}
\usage{
	network.AIC(network,Pm=NULL,k=length(nodes(network$graph)),verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{a nem object (e.g. 'pairwise')}
  \item{Pm}{prior over models (n x n matrix). If NULL, then a matrix of 0s is assumed}
  \item{k}{penalty per parameter in the AIC/BIC calculation. k = 2 for classical AIC}
  \item{verbose}{print out the result}  
}

\details{
	For k = log(n) the BIC (Schwarz criterion) is computed. Usually this function is not called directly but from \code{nemModelSelection}
}

\value{
	AIC/BIC value
}

\author{Holger Froehlich}


\seealso{\code{\link{nemModelSelection}}}
\examples{
   data("BoutrosRNAi2002") 
   D = BoutrosRNAiDiscrete[,9:16]
   control = set.default.parameters(unique(colnames(D)), para=c(0.13,0.05))
   res1 <- nem(D, control=control)
   network.AIC(res1)
   control$lambda=100 # enforce sparsity
   res2 <- nem(D,control=control)
   network.AIC(res2)
}
\keyword{graphs}
