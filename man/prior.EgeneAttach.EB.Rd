\name{prior.EgeneAttach.EB}
\alias{prior.EgeneAttach.EB}

\title{Initialize E-gene attachment prior for empirical Bayes}
\description{This function computes an initial prior for the attachment of E-genes based on observed log-pvalue densities or log-odds ratios. This prior can then be updated within an empirical Bayes procedure (MC.EMiNEM).}
\usage{
prior.EgeneAttach.EB(ratioMat)
}

\arguments{
  \item{ratioMat}{data matrix with experiments in the columns (log-odds ratios or log-pvalue densities)}  
}
\value{
  |E-genes| x (|S-genes| + 1) matrix with prior E-gene attachment probabilities. The last column denotes the virtual 'null' S-gene, which is there to filter E-genes that have no obvious attachment to any of the real S-genes. 
}

\author{Theresa Niederberger, Holger Froehlich}

\references{
Niederberger, T.; Etzold, S.; Lidschreiber, M; Maier, K.; Martin, D.; Fr\"ohlich, H.; Cramer, P.; Tresch, A., MC Eminem Maps the Interaction Landscape of the Mediator, PLoS Comp. Biol., 2012, submitted.
}
\seealso{\code{\link{set.default.parameters}}, \code{\link{nem}}}
\examples{
# only for test purposes
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDens
   control = set.default.parameters(unique(colnames(D)), Pe=prior.EgeneAttach.EB(D), mcmc.nsamples=100, mcmc.nburnin=50, type="CONTmLLBayes")  # these are *not* realistic values
   res <- nem(D,inference="mc.eminem", control=control)        
   if(require(Rgraphviz))    
    plot(res)   
}
\keyword{graphs}
\keyword{models}
