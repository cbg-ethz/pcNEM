\name{getDensityMatrix}
\alias{getDensityMatrix}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate density matrix from raw p-value matrix}
\description{
  Fit a 3 component BUM model to each column of a raw p-value matrix.
}
\usage{
	getDensityMatrix(Porig, dirname=NULL, startab=c(0.3,10), startlam=c(0.6,0.1,0.3), tol=1e-4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Porig}{matrix of raw p-values}
\item{dirname}{name of a directory to save histograms and QQ-plots to. If dirname=NULL, then the plots are made to the screen, and after each fit the user is asked to press a key in order to continue.}
\item{startab}{start values for alpha and beta parameter}
\item{startlam}{start values for mixing coefficients}
\item{tol}{convergence tolerance: If the absolute likelihood ratio -1 becomes smaller than this value, then the EM algorithm is supposed to be converged.}
}
\details{
	The BUM density model consists of 3 components: \eqn{f(x) = lambda_1 + lambda_2*dbeta(x,alpha,1) + lambda_3*dbeta(x,1,beta)}. The mixing coefficients and the parameters alpha and beta are fitted together via an EM algorithm.
}
\value{
 	log-density matrix of same dimensions as Porig: The log-densities can be interpreted as log signal-to-noise ratios. A value > 0 means higher signal than noise, and a value < 0 a higher noise than signal.
}
\note{Note the difference to the previous package version: the LOG-density is returned now!}
\author{ Holger Froehlich }
\keyword{models}