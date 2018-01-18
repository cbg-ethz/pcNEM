\name{NiederbergerMediator2012}
\alias{NiederbergerMediatorLods}
\alias{NiederbergerMediatorLogFC}
\alias{NiederbergerMediatorPVals}
\alias{NiederbergerMediator2012}
\docType{data}
\title{Expression measurements upon perturbation of Mediator subunits}
\description{
This dataset combines expression data from 4 gene perturbation studies in S. cerevisae, including the perturbations of the 9 Mediator 
subunits Med2, Med7 (N-terminus), Med7 (C-terminus), Med10, Med15, Med19, Med20, Med21 and Med31. Med10 and Med21 have been combined due to indistinguishable results.
}
\usage{data(NiederbergerMediator2012)}
\format{
  NiederbergerMediatorLods:	data matrix(log-odds ratios):		2578 x 8
  NiederbergerMediatorLogFC: 	data matrix(log2 fold-changes):		2578 x 8
  NiederbergerMediatorPVals: 	data matrix(p-values):			2578 x 8
}
\details{
The data has been preprocessed using the R/Bioconductor package limma.
The results have been filtered as follows: The original data contains a double-knockout Med7N/Med31, which can't be assessed by EMiNEM and has thus been removed. 
Analysis of the data revealed genes showing batch-specific effects, they have been excluded from further analysis. Finally, genes that do not react to any 
perturbation have been removed since they do not contain additional information. For a more detailed description of the individual studies and datasets as well as 
further information on data processing and results, please refer to Niederberger et al, 2012. The raw data is available on Array Express (accession number [provided upon publication]).

\code{NiederbergerMediatorLods}:	the log-odds ratios
\code{NiederbergerMediatorLogFC}:	the log2 fold-changes
\code{NiederbergerMediatorPVals}:	the corresponding p-values

}
\references{
Niederberger T., Etzold S., Lidschreiber M., Maier K.C., Martin D.E., Fr\"o hlich H., Cramer P., Tresch A., MC EMiNEM Maps the Interaction Landscape of the Mediator. PLoS Computational Biology, 8(6): e1002568, 2012-
}
\examples{
    data(NiederbergerMediator2012)
    dim(NiederbergerMediatorLods)
    dim(NiederbergerMediatorLogFC)
    dim(NiederbergerMediatorPVals)
}
\keyword{datasets}

