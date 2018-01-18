\name{quicknem}
\alias{quicknem}

\title{Quick run of Nested Effects Models inference}
\description{Interface to learn NEM models from data}
\usage{
quicknem(D,type="CONTmLLDens",inference="nem.greedy",controls.name=NULL,contrasts=NULL,normalize=FALSE,cutoff=0.05,DIR="bum",plot=TRUE,bootstrap=0,...) 
}

\arguments{
  \item{D}{ExpressionSet object or data matrix with raw or normalized expression data.}
  \item{type}{Parameter estimation, either \code{mLL}, \code{FULLmLL}, \code{CONTmLL}, \code{CONTmLLBayes}, \code{CONTmLLMAP}, \code{depn}.}
  \item{inference}{\code{search} to use exhaustive enumeration, \code{triples} for triple-based inference, \code{pairwise} for the pairwise heuristic, \code{ModuleNetwork} for the module based inference, \code{nem.greedy} for greedy hillclimbing, \code{nem.greedyMAP} for alternating MAP optimization using log odds or log p-value densities}  
  \item{controls.name}{Pattern to search for in the columnnames of D. Defines which columns in D should be regarded as controls.}
  \item{contrasts}{String defining the contrasts to estimate via limma}
  \item{normalize}{boolean value, should quantile normalization be performed}
  \item{cutoff}{P-value cutoff for differential expression using adjusted p.values from limma.}
  \item{DIR}{Directory name, where additional informative plots should be stored. Created if not present.}
  \item{plot}{Should the inferred network be plotted?}
  \item{bootstrap}{Integer defining the number of bootstrapping samples to be performed. Defaults to 0.}
  \item{...}{other arguments to pass}
}
\details{
	Wrapper function for call of nem inference. Extracts differential genes for given contrasts and infers a NEM - graph for the given inference type.
\describe{
	\item{D}{Is either an ExpressionSet Object or a matrix/data.frame containing the expression values from the siRNA knockdown experiments. If an ExpressionSet, the data is extracted via exprs(ExpressionSet). The knockdowns must be in the columns, the measured effect genes in the rows of the expression matrix.}
	\item{type}{\code{mLL} or \code{FULLmLL} or \code{CONTmLL} or \code{CONTmLLBayes} or \code{CONTmLLMAP} or \code{depn}. \code{CONTmLLDens} and \code{CONTmLLRatio} are identical to \code{CONTmLLBayes} and \code{CONTmLLMAP} and are still supported for compatibility reasons. \code{mLL} and \code{FULLmLL} are used for binary data (see \code{BoutrosRNAiDiscrete}) and \code{CONTmLL} for a matrix of effect probabilities. \code{CONTmLLBayes} and \code{CONTmLLMAP} are used, if log-odds ratios, p-value densities or any other model specifies effect likelihoods. \code{CONTmLLBayes} refers to an inference scheme, were the linking positions of effect reporters to network nodes are integrated out, and \code{CONTmLLMAP} to an inference scheme, were a MAP estimate for the linking positions is calculated. \code{depn} indicates Deterministic Effects Propagation Networks (DEPNs).}
	\item{inference}{Type of network reconstruction. \code{search} enumerates all possible networks. Set to \code{triples}, \code{pairwise}, \code{ModuleNetwork}, \code{nem.greedy} or \code{nem.greedyMAP} for heuristic search of the network.}
	\item{controls.name}{Defines a pattern to search for in the column names of D, which describes the control experiments. Each remaining experiment is then compared via limma to these controls by defining the appropriate contrasts. If NULL, then \code{controls.name} must be given, except for using \code{type="depn"}, where neither \code{controls.name} nor \code{contrasts} needs to be defined.}
	\item{contrasts}{Defines the contrasts of interest that should be used for the limma analysis. If NULL, then \code{controls.name} must be given, except for using \code{type="depn"}, where neither \code{controls.name} nor \code{contrasts} needs to be defined.}
	\item{DIR}{In case of \code{type="CONTmLLDens"} or \code{type="CONTmLLBayes"} some additional plots for the BUM model fits are created and stored here.}
}
}
\value{
  \item{graph}{the inferred directed graph (graphNEL object)}   
  \item{mLL}{log posterior marginal likelihood of final model} 
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions}   
  \item{selected}{selected E-gene subset}
  \item{LLperGene}{likelihood per selected E-gene}
  \item{control}{hyperparameter as in function call}
  \item{bootstrap}{Integer number defining how many bootstrap samples should be drawn. If 0, no bootstrapping will be performed. Else, \code{\link{nem.bootstrap}} will be called internally.} 
}

\author{Christian Bender, Holger Froehlich, Florian Markowetz}


\seealso{\code{\link{nem}}, \code{\link{set.default.parameters}}, \code{\link{nemModelSelection}}, \code{\link{nem.jackknife}}, \code{\link{nem.bootstrap}}, \code{\link{nem.consensus}}, \code{\link{local.model.prior}}, \code{\link{plot.nem}}}
\examples{
\dontrun{
data(BoutrosRNAi2002)
exps <- colnames(BoutrosRNAiExpression)
res <- quicknem(BoutrosRNAiExpression,controls="control")
res <- quicknem(BoutrosRNAiExpression,controls="control",type="CONTmLLRatio")
res <- quicknem(BoutrosRNAiExpression,controls="control",type="CONTmLLRatio",inference="ModuleNetwork")
contrasts <- c("rel-control","rel-LPS","key-control","key-LPS","tak-control","tak-LPS","mkk4hep-control","mkk4hep-LPS")
res <- quicknem(BoutrosRNAiExpression,contrasts=contrasts)

data(SahinRNAi2008)
dat <- dat.unnormalized #[,sample(1:17,5)]
res <- quicknem(dat,type="depn")
}
}
\keyword{graphs}
\keyword{models}
