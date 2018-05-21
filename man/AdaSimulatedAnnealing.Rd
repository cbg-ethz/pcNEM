\name{AdaSimulatedAnnealing}
\alias{AdaSimulatedAnnealing}
\title{Adaptive simulated annealing based inference algorithm}
\description{Internal function for pc-NEM for inference of graph structure and noise paramters (optional) using adaptive simulated annealing. The code for adaptive simulated annealing has been developed using the code from the paper titled 'Partition MCMC for Inference on Acyclic Digraphs' by Jack Kuipers & Giusi Moffa which is further based on the code from the Dortmund course programmed by Miriam Lohr. The code has been modified to match it to the NEMs framework. The scoring function has been changed to compute the likelihood according to our model. Further, the update steps have been modified.
}
\usage{
nem.AdaSA(n,phi,D,control)
}
\arguments{
  \item{n}{Number of S-genes}
  \item{phi}{The starting S-gene model}
  \item{\code{D}}{The data set. The rows refer to the E-genes and the columns correspond to the knockdown experiment. \code{D} has to be binary and the column names should match the row names of the perturbation map}
  \item{control}{Control parameters set using set.default.parameters. Object stores the information on number of iterations, temperature, adaptation rate, ideal acceptance rate, perturbation map and option to infer noise parameters.}
}
\details{
  The algorithm stochastically alternates between two distinct spaces, one for the graphs and one for the error rates.
}
\value{
\item{graphs}{Graph chain}
\item{LLscores}{Log-likelihood score chain}
\item{Temp}{Temperature chain}
\item{AcceptRate}{Acceptance rate chain}
\item{AlphaVals}{False positive rates chain}
\item{BetaVals}{False negative rates chain}
\item{maxLLscore}{Log-likelihood value of best graph}
\item{maxDAG}{DAG with the maximum log-likelihood value}
\item{typeI}{MLE of false positive rate}
\item{typeII}{MLE of false negative rate}
\item{minIter}{Minimum steps taken to reach the maximum log-likeliood value}
\item{transformAR}{Log transformed acceptance rate}
}
\references{Srivatsa S, Kuipers J, Schmich F, Eicher S, Emmenlauer M, Dehio C, Beerenwinkel N, Improved pathway reconstruction from RNA interference screens by exploiting off-target effects, ISMB, 2018}
\author{Sumana Srivatsa}
\note{Can be modified for continuous data if needed.}

\keyword{graphs}
\keyword{models}
