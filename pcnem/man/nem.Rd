\name{nem}
\alias{nem}
\alias{print.nem}
\alias{print.nem.greedy}
\alias{print.nem.greedyMAP}
\alias{print.pairwise}
\alias{print.triples}
\alias{print.ModuleNetwork}
\alias{print.score}
\alias{print.nem.BN}
\alias{print.mc.eminem}
\alias{print.dynoNEM}

\title{Nested Effects Models - main function}
\description{The main function to perform model learning from data}
\usage{
nem(D,inference="nem.greedy",models=NULL,control=set.default.parameters(setdiff(unique(colnames(D)),"time")), verbose=FALSE)

\method{print}{nem}(x, ...)
}

\arguments{
  \item{D}{data matrix with experiments in the columns (binary or continious)}
  \item{inference}{\code{search} to use exhaustive enumeration, \code{triples} for triple-based inference, \code{pairwise} for the pairwise heuristic, 
  \code{ModuleNetwork} for the module based inference proposed in Fr\"ohlich et al. 2008, \code{ModuleNetwork.orig} for the module based inference proposed in Fr\"ohlich et al. 2007, \code{nem.greedy} for greedy hillclimbing, \code{nem.greedyMAP} for alternating MAP optimization using log odds or log p-value densities, \code{mc.eminem} for EM based inference using log odds or log p-value densities, \code{BN.greedy}, \code{BN.exhaustive} for a conventional Bayesian Network treatment using binomial or normal distribution assumptions, \code{dynoNEM} for MCMC based inference from time series data, \code{mc.eminem} for EM based inference}
  \item{models}{a list of adjacency matrices for model search. If NULL, an  exhaustive enumeration of all possible models is performed.}
  \item{control}{list of parameters: see \code{set.default.parameters}}
  \item{verbose}{do you want to see progression statements? Default: TRUE}

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\details{  
   If parameter Pm != NULL and parameter lambda == 0, a Bayesian approach to include prior knowledge is used. Alternatively, the regularization parameter lambda can be tuned in a model selection step via the function \code{nemModelSelection} using the BIC criterion. 
   If automated subset selection of effect reporters is used (default), the regularization parameter delta can be tuned via the BIC model selection criterion. Per default it is fixed to 1 / (no. S-genes + 1). 

  The function \code{plot.nem} plots the inferred phenotypic hierarchy as a directed graph, the likelihood distribution of the models (only for exhaustive search) or the posterior position of the effected genes.  
}
\value{
  \item{graph}{inferred directed S-gene graph (graphNEL object)}   
  \item{mLL}{log posterior marginal likelihood of model(s)} 
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions}   
  \item{selected}{selected E-gene subset}
  \item{LLperGene}{likelihood per selected E-gene}
  \item{avg}{in case of MCMC: posterior mean S-gene graph (edge weighted adjacency matrix)}
  \item{control}{hyperparameter as in function call}
  
  For inference = "mc.eminem" the following additional values are returned:
  \item{local.maxima}{local maxima of the EM procedure}  
  \item{graphs.sampled}{sampled graphs}
  \item{EB}{samples of the empirical Bayes prior}
  \item{acc_list}{list that indicates whether the corresponding sampled S-gene graph has been accepted (new local maximum (1), same local maximum (0)) or rejected(-1) in the MCMC sampling process - length(acc_list)=mcmc.nsamples + mcmc.nburnin}
}
\references{

Markowetz, F.; Bloch, J. & Spang, R., Non-transcriptional Pathway Features Reconstructed from Secondary Effects of RNA interference. Bioinformatics, 2005, 21, 4026 - 4032\\

Markowetz, F.; Kostka, D.; Troyanskaya, O. & Spang, R., Nested Effects Models for High-dimensional Phenotyping Screens. Bioinformatics, 2007, 23, i305 - i312\\

Fr\"ohlich, H.; Fellmann, M.; S\"ultmann, H.; Poustka, A. & Beissbarth, T. Large Scale Statistical Inference of Signaling Pathways from RNAi and Microarray Data. BMC Bioinformatics, 2007, 8, 386\\

Fr\"ohlich, H.; Fellmann, M.; S\"ultmann, H.; Poustka, A. & Beissbarth, T. Estimating Large Scale Signaling Networks through Nested Effect Models with Intervention Effects from Microarray Data. Bioinformatics, 2008, 24, 2650-2656\\ 

Tresch, A. & Markowetz, F., Structure Learning in Nested Effects Models Statistical Applications in Genetics and Molecular Biology, 2008, 7\\

Zeller, C.; Fr\"ohlich, H. & Tresch, A., A Bayesian Network View on Nested Effects Models EURASIP Journal on Bioinformatics and Systems Biology, 2009, 195272\\

Fr\"ohlich, H.; Tresch, A. & Beissbarth, T., Nested Effects Models for Learning Signaling Networks from Perturbation Data. Biometrical Journal, 2009, 2, 304 - 323\\

Fr\"ohlich, H.; Sahin, \"O.; Arlt, D.; Bender, C. & Beissbarth, T. Deterministic Effects Propagation Networks for Reconstructing Protein Signaling Networks from Multiple Interventions. BMC Bioinformatics, 2009, 10, 322\\

Fr\"ohlich, H.; Praveen, P. & Tresch, A., Fast and Efficient Dynamic Nested Effects Models. Bioinformatics, 2011, 27, 238-244\\

Niederberger, T.; Etzold, S.; Lidschreiber, M; Maier, K.; Martin, D.; Fr\"ohlich, H.; Cramer, P.; Tresch, A., MC Eminem Maps the Interaction Landscape of the Mediator, PLoS Comp. Biol., 8(6): e1002568, 2012.

}

\author{Holger Froehlich, Florian Markowetz}


\seealso{\code{\link{set.default.parameters}}, \code{\link{nemModelSelection}}, \code{\link{nem.jackknife}}, \code{\link{nem.bootstrap}}, \code{\link{nem.consensus}}, \code{\link{local.model.prior}}, \code{\link{plot.nem}}}
\examples{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   control = set.default.parameters(unique(colnames(D)), para=c(0.13, 0.05))   
   res1 <- nem(D,inference="search", control=control)
   res2 <- nem(D,inference="pairwise", control=control)
   res3 <- nem(D,inference="triples", control=control)
   res4 <- nem(D,inference="ModuleNetwork", control=control)
   res5 <- nem(D,inference="nem.greedy", control=control)        
   res6 = nem(BoutrosRNAiLods, inference="nem.greedyMAP", control=control)
   

   par(mfrow=c(2,3))
   plot.nem(res1,main="exhaustive search")
   plot.nem(res2,main="pairs")
   plot.nem(res3,main="triples")
   plot.nem(res4,main="module network")
   plot.nem(res5,main="greedy hillclimber")      
   plot.nem(res6,main="alternating MAP optimization")
}
\keyword{graphs}
\keyword{models}
