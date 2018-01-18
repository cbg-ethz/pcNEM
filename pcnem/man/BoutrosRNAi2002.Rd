\name{BoutrosRNAi2002}
\alias{BoutrosRNAiExpression}
\alias{BoutrosRNAiDiscrete}
\alias{BoutrosRNAiDens}
\alias{BoutrosRNAiLods}
\alias{BoutrosRNAiLogFC}
\alias{BoutrosRNAi2002}
\docType{data}
\title{RNAi data on Drosophila innate immune response}
\description{
Data from a study on innate immune response in \emph{Drosophila} (Boutros
et al, 2002). Selectively removing signaling components by RNAi blocked induction of all, or only parts, of the
transcriptional response to LPS. The nested structure of perturbation effects allows to reconstruct a
branching in the Imd pathway.
}
\usage{data(BoutrosRNAi2002)}
\format{
  BoutrosRNAiExpression: data matrix: 14010 x 16\\
  BoutrosRNAiDiscrete:   binary matrix:  68 x 16\\
  BoutrosRNAiDens:		data matrix: 68 x 4\\
  BoutrosRNAiLods: 	    data matrix: 68 x 4\\
  BoutrosRNAiLogFC:		data matrix: 68 x 4
}
\details{
The dataset consists of 16 Affymetrix-microarrays: 4 replicates of control
experiments without LPS and without RNAi (negative controls), 4 replicates of
expression profiling after stimulation with LPS but without RNAi (positive
controls), and 2 replicates each of expression profiling after applying LPS and
silencing one of the four candidate genes tak, key, rel, and mkk4/hep. 

\code{BoutrosRNAiExpression}: For preprocessing we performed normalization on probe level using a variance 
stabilizing transformation (Huber et al, 2002), and probe set summarization
using a median polish fit of an additive model (Irizarry et al, 2003).

\code{BoutrosRNAiDiscrete}: contains only the 68 genes more than two-fold up-regulated between negative and positive controls.
The continuous expression values are discretized to \code{1} (effect: closer to negative controls) and \code{0} (no effect: closer to positive controls).

\code{BoutrosRNAiDens}: log \eqn{p}-value density matrix for the 68 genes with more than two-fold up-regulated between negative and positive controls.

\code{BoutrosRNAiLods}: B-value matrix for the 68 genes with more than two-fold up-regulated between negative and positive controls.

\code{BoutrosRNAiLogFC}: matrix with log fold changes

}
\references{
  Boutros M, Agaisse H, Perrimon N, Sequential activation of signaling pathways
  during innate immune responses in Drosophila. Developmental
  Cell. 3(5):711-722, 2002
}
\seealso{\code{\link{nem.discretize}}}
\examples{
    data("BoutrosRNAi2002")
    dim(BoutrosRNAiExpression)
    dim(BoutrosRNAiDiscrete)      
}
\keyword{datasets}
