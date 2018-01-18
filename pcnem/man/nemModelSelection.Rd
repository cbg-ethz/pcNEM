\name{nemModelSelection}
\alias{nemModelSelection}

\title{Model selection for nested effect models}
\description{Infers models with different regularization constants, compares them via the BIC or AIC criterion and returns the highest scoring one}
\usage{
nemModelSelection(lambdas,D,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))),verbose=TRUE,...)
}

\arguments{
  \item{lambdas}{vector of regularization constants}
  \item{D}{data matrix with experiments in the columns (binary or continious)}
  \item{inference}{\code{search} to use exhaustive enumeration, \code{triples} for triple-based inference, \code{pairwise} for the pairwise heuristic, \code{ModuleNetwork} for the module based inference, \code{nem.greedy} for greedy hillclimbing, \code{nem.greedyMAP} for alternating MAP optimization using log odds or log p-value densities}
  \item{models}{a list of adjacency matrices for model search. If NULL, an exhaustive enumeration of all possible models is performed.}
  \item{control}{list of parameters: see \code{set.default.parameters}}
  \item{verbose}{do you want to see progression statements? Default: TRUE}
  \item{...}{other arguments to pass to function \code{nem} or \code{network.AIC}}
}
\details{
  \code{nemModelSelection} internally calls \code{nem} to infer a model with a given regularization constant. The comparison between models is based on the BIC or AIC criterion, depending on the parameters passed to \code{network.AIC}.
  
}
\value{
  nem object
}

\author{Holger Froehlich}


\seealso{\code{\link{set.default.parameters}}, \code{\link{nem}}, \code{\link{network.AIC}}}
\examples{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]   
   hyper = set.default.parameters(unique(colnames(D)), para=c(0.13, 0.05), Pm=diag(4))
   res <- nemModelSelection(c(0.1,1,10), D, control=hyper)      
   
   plot.nem(res,main="highest scoring model")      
}
\keyword{graphs}
\keyword{models}
