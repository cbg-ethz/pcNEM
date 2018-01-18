\name{nem.jackknife}
\alias{nem.jackknife}
\alias{print.nem.jackknife}

\title{Jackknife for nested effect models}
\description{   
 Assesses the statistical stability of a network via a jackknife procedure: Each S-gene is left out once and the network reconstructed on the remaining ones. The relative frequency of each edge to appear in n-1 jackknife samples is returned.
}
\usage{
nem.jackknife(D, thresh=0.5, inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))), verbose=TRUE)

\method{print}{nem.jackknife}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix with experiments in the columns (binary or continious)}
  \item{thresh}{only edges appearing with a higher frequency than "thresh" are returned}
  \item{inference}{\code{search} to use exhaustive enumeration, \code{triples} for triple-based inference, \code{pairwise} for the pairwise heuristic, \code{ModuleNetwork} for the module based inference, \code{nem.greedy} for greedy hillclimbing, \code{nem.greedyMAP} for alternating MAP optimization using log odds or log p-value densities}
  \item{models}{a list of adjacency matrices for model search. If NULL, an  exhaustive enumeration of all possible models is performed.}
  \item{control}{list of parameters: see \code{set.default.parameters}}
  \item{verbose}{do you want to see progression statements? Default: TRUE}

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\details{
Calls \code{\link{nem}} or \code{\link{nemModelSelection}} internally, depending on whether or not parameter lambda is a vector and parameter Pm != NULL.
}
\value{
nem object with edge weights being the jackknife probabilities
}

\author{Holger Froehlich}

\seealso{\code{\link{nem.bootstrap}}, \code{\link{nem.consensus}}, \code{\link{nem}}, \code{\link{nemModelSelection}}}
\examples{
\dontrun{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   nem.jackknife(D, control=set.default.parameters(unique(colnames(D)), para=c(0.13,0.05)))         
}
}

\keyword{models}
