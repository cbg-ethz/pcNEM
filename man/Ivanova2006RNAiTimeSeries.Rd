\name{Ivanova2006RNAiTimeSeries}
\alias{Ivanova2006RNAiTimeSeries}
\alias{dat}
\docType{data}
\title{Perturbation Time Series}
\description{
The dataset consists of time series microarray measurements (only 1 replicate, 6 time points) for RNAi knock-downs of each of the six genes Nanog, Oct4, Sox2, Esrrb, Tbx3, and Tcl1 in mouse embryonic stem cells.
}
\usage{data(Ivanova2006RNAiTimeSeries)}
\format{
 dat: 8 x 122 x 6 array of discretized time series data for E-genes showing a abs(logFC) > 1 in at least one time point (see Anchang et al., 2009).
}
\details{
Mouse embryonic stem cells
(ESCs) were grown in the presence of leukemia inhibitory factor
(LIF), thus retaining their undifferentiated self-renewing state (positive
controls). Differentiation associated changes in gene expression
were measured by replacing LIF with retinoic acid (RA), thus inducing
differentiation of stem cells (negative controls). RNAi silencing
of the 6 afore mentioned genes was done in (LIF+, RA-) cell culturs
to investigate their potential for induced cell differentiation. Microarray
expression measurements at 6 - 7 time points in one-day
intervals were taken for the two controls (positive and negative) and
the six RNAi assays.
The dataset by Ivanova et al. (2006) was measured once on
Affymetrix MOE430A and once on MOE430B chips. The authors
normalized both datasets via the MAS 5.0 software (Hubbell et al., 2002). Here, only the Affymetrix MOE430A chip series is included.
}
\references{
Ivanova, N., Dobrin, R., Lu, R., Kotenko, I., Levorse, J., DeCoste, C., Schafer, X., Lun, Y., and Lemischka, I. R. (2006). Dissecting self-renewal in stem cells with RNA interference. Nature, 442(7102):533-538.\\

Anchang, B., Sadeh, M. J., Jacob, J., Tresch, A., Vlad, M. O., Oefner, P. J., and Spang, R. (2009). Modeling the temporal interplay of molecular signaling and gene expression by using dynamic nested effects models. Proceedings of the National Academy of Sciences of the United States of America, 106(16):6447-52.
}
\seealso{\code{\link{Ivanova2006RNAiTimeSeries}}}
\examples{
    data("Ivanova2006RNAiTimeSeries")
    dim(dat)          
}
\keyword{datasets}
