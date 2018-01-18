\name{sim.intervention}
\alias{sim.intervention}
\alias{sim.intervention}

\title{Simulate interventions in a Nested Effects Model}
\description{   
 Simulates a knock-down of a list of network nodes and returns the network nodes and effect reporters, where effects are expected.
}
\usage{
sim.intervention(x, int, T=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{nem object}
  \item{int}{a character vector of nodes in the network}
  \item{T}{number of time steps to be simulated (only for dynoNEMs)}  
}
\value{
list with two slots:
\item{Sgenes.effected}{effected network nodes (by intervention)}
\item{Egenes.effected}{most likely effected downstream effect reporters (by intervention)}
\item{Egenes.effected.prob}{probability of each reporter to be effected by the respective intervention experiment.}
}

\author{Holger Froehlich}

\examples{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   res = nem(D, control=set.default.parameters(unique(colnames(D)), para=c(0.13,0.05))) 
   sim.intervention(res, "rel") # simulate knock-down of rel
}

\keyword{models}
