\name{nem.cont.preprocess}
\alias{nem.cont.preprocess}
\title{Calculate classification probabilities of perturbation data according to control experiments}
\description{Calculates probabilities of data to define effects of interventions with respect to wildtype/control measurements}
\usage{
nem.cont.preprocess(D,neg.control=NULL,pos.control=NULL,nfold=2, influencefactor=NULL, empPval=.05, verbose=TRUE)
}
\arguments{
  \item{D}{matrix with experiments as columns and effect reporters as rows}
  \item{neg.control}{either indices of columns in \code{D} or a matrix with the same number of rows as \code{D}}
  \item{pos.control}{either indices of columns in \code{D} or a matrix with the same number of rows as \code{D}}
  \item{nfold}{fold-change between neg. and pos. controls for selecting effect reporters. Default: 2}
  \item{influencefactor}{factor multiplied onto the probabilities, so that all negative control genes are treated as influenced, usually automatically determined}
  \item{empPval}{empirical p-value cutoff for effects if only one control is available}
  \item{verbose}{Default: TRUE}
}
\details{
  Determines the empirical distributions of the controls and calculates the probabilities of pertubartion data to belong to the control distribution(s).
}
\value{
\item{dat}{data matrix}
\item{pos}{positive controls [in the two-controls setting]}
\item{neg}{negative controls [in the two-controls setting]}
\item{sel}{effect reporters selected [in the two-controls setting]}
\item{prob.influenced}{probability of a reporter to be influenced}
\item{influencefactor}{factor multiplied onto the probabilities, so that all negative control genes are treated as influenced}
}
\references{Markowetz F, Bloch J, Spang R, Non-transcriptional pathway features reconstructed from secondary effects of RNA interference, Bioinformatics, 2005}
\author{Florian Markowetz}
\note{preliminary! will be developed to be more generally applicable}

\seealso{\code{\link{BoutrosRNAi2002}}}
\examples{
   data("BoutrosRNAi2002")
   preprocessed <- nem.cont.preprocess(BoutrosRNAiExpression,neg.control=1:4,pos.control=5:8)
}
\keyword{graphs}
\keyword{models}
