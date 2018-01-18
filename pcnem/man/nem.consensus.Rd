\name{nem.consensus}
\alias{nem.consensus}
\alias{print.nem.consensus}

\title{Statistically stabile nested effects models}
\description{Performs bootstrapping (resampling with replacement) on E-genes and jackknife on S-genes to assess the statistical stability of networks. Only edges appearing with a higher frequency than a predescribed threshold in both procedures are regarded as statistical stable and appear in the so-called consensus network.}
\usage{
nem.consensus(D,thresh=0.5, nboot=1000,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))), verbose=TRUE)

\method{print}{nem.consensus}(x, ...)
}

\arguments{  
  \item{D}{data matrix with experiments in the columns (binary or continous)}
  \item{thresh}{only edges appearing with a higher frequency than "thresh" in both, bootstrap and jackknife procedure, are regarded as statistically stable and trust worthy}
  \item{nboot}{number of bootstrap samples desired}
  \item{inference}{\code{search} to use exhaustive enumeration, \code{triples} for triple-based inference, \code{pairwise} for the pairwise heuristic, \code{ModuleNetwork} for the module based inference, \code{nem.greedy} for greedy hillclimbing, \code{nem.greedyMAP} for alternating MAP optimization using log odds or log p-value densities}
  \item{models}{a list of adjacency matrices for model search. If NULL, an  exhaustive enumeration of all possible models is performed.}
  \item{control}{list of parameters: see \code{set.default.parameters}}
  \item{verbose}{do you want to see progression statements? Default: TRUE}

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\details{
  Calls \code{\link{nem}} or \code{\link{nemModelSelection}} internally, depending on whether or not lambda is a vector and Pm != NULL.  
}
\value{
	consensus network (nem object)
}

\author{Holger Froehlich}


\seealso{\code{\link{nem.bootstrap}}, \code{\link{nem.jackknife}}, \code{\link{nem.calcSignificance}}, \code{\link{nem}}}
\examples{
\dontrun{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]   
   nem.consensus(D, control=set.default.parameters(unique(colnames(D)), para=c(0.13,0.05)))            
}
}
\keyword{graphs}
\keyword{models}
