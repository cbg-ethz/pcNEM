\name{nem.discretize}
\alias{nem.discretize}
\title{Discretize perturbation data according to control experiments}
\description{discretizes raw data to define effects of interventions with respect to wildtype/control measurements}
\usage{
nem.discretize(D,neg.control=NULL,pos.control=NULL,nfold=2,cutoff=0:10/10, pCounts=20, empPval=.05, verbose=TRUE)
}
\arguments{
  \item{D}{matrix with experiments as columns and effect reporters as rows}
  \item{neg.control}{either indices of columns in \code{D} or a matrix with the same number of rows as \code{D}}
  \item{pos.control}{either indices of columns in \code{D} or a matrix with the same number of rows as \code{D}}
  \item{nfold}{fold-change between neg. and pos. controls for selecting effect reporters. Default: 2}
  \item{cutoff}{a (vector of) cutoff value(s) weighting the pos. controls versus the neg. controls. Default: 0:10/10} 
  \item{pCounts}{pseudo-counts to guard against unreasonable low error estimates}
  \item{empPval}{empirical p-value cutoff for effects if only one control is available}
  \item{verbose}{Default: TRUE}
}
\details{
  Chooses cutoff such that separation between negative and positive controls becomes optimal.
}
\value{
\item{dat}{discretized data matrix}
\item{pos}{discretized positive controls [in the two-controls setting]}
\item{neg}{discretized negative controls [in the two-controls setting]}
\item{sel}{effect reporters selected [in the two-controls setting]}
\item{cutoff}{error rates for different cutoff values [in the two-controls setting]}
\item{para}{estimated error rates [in the two-controls setting]}
}
\references{Markowetz F, Bloch J, Spang R, Non-transcriptional pathway features reconstructed from secondary effects of RNA interference, Bioinformatics, 2005}

\author{Florian Markowetz}
\note{preliminary! will be developed to be more generally applicable}

\seealso{\code{\link{BoutrosRNAi2002}}}
\examples{
   # discretize Boutros data as in
   # Markowetz et al, 2005
   data("BoutrosRNAi2002")
   disc <- nem.discretize(BoutrosRNAiExpression,neg.control=1:4,pos.control=5:8,cutoff=.7)
   stopifnot(disc$dat==BoutrosRNAiDiscrete[,9:16])   
}
\keyword{graphs}
\keyword{models}
