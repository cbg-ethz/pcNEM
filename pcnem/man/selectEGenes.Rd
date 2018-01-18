\name{getRelevantEGenes}
\alias{getRelevantEGenes}
\alias{selectEGenes}
\alias{filterEGenes}

\title{Automatic selection of most relevant effect reporters}
\description{   
   1. A-priori filtering of effect reporters/E-genes: Select effect reporters, which show a pattern of differential expression across experiments that is expected to be non-random. 
   2. Automated effect reporters subset selection: Select those effect reporters, which have the highest likelihood under the given network hypothesis.
}
\usage{
filterEGenes(Porig, D, Padj=NULL, ntop=100, fpr=0.05, adjmethod="bonferroni", cutoff=0.05)

getRelevantEGenes(Phi, D, control, nEgenes=min(10*nrow(Phi), nrow(D)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
For method filterEGenes:
  \item{Porig}{matrix of raw p-values, typically from the complete array}
  \item{D}{ data matrix. Columns correspond to the nodes in the silencing scheme. Rows are effect reporters. }
  \item{Padj}{matrix of false positive rates. If not, provided Benjamini-Hochbergs method for false positive rate computation is used.}
  \item{ntop}{number of top genes to consider from each knock-down experiment}
  \item{fpr}{significance cutoff for the FDR}
  \item{adjmethod}{adjustment method for pattern p-values}
  \item{cutoff}{significance cutoff for patterns}
  
For method getRelevantEGenes:
  \item{Phi}{ adjacency matrix with unit main diagonal }    
  \item{control}{list of parameters: see \code{set.default.parameters}}
  \item{nEgenes}{ no. of E-genes to select}  
}
\details{
	The method filterEGenes performs an a-priori filtering of the complete microarray. It determines how often E-genes are expected to be differentially expressed across experiments just randomly. According to this only E-genes are chosen, which show a pattern of differential expression more often than can be expected by chance.

	The method getRelevantEGenes looks for the E-genes, which have the highest likelihood under the given network hypothesis. In case of the scoring type "CONTmLLBayes" these are all E-genes which have a positive contribution to the total log-likelihood. In case of type "CONTmLLMAP" all E-genes not assigned to the "null" S-gene are returned. This involves the prior probability delta/no. S-genes for leaving out an E-gene.
	For all other cases ("CONTmLL", "FULLmLL", "mLL") the nEgenes E-genes with the highest likelihood under the given network hypothesis are returned.
}
\value{
  \item{I}{index of selected E-genes}
  \item{dat}{subset of original data according to I }
  \item{patterns}{significant patterns}
  \item{nobserved}{no. of cases per observed pattern}  



  \item{selected}{selected E-genes}
  \item{mLL}{marginal likelihood of a phenotypic hierarchy}
  \item{pos}{posterior distribution of effect positions in the hierarchy}
  \item{mappos}{Maximum a posteriori estimate of effect positions}
  \item{LLperGene}{likelihood per selected E-gene}
}

\author{Holger Froehlich}

\seealso{\code{\link{nem}}, \code{\link{score}}, \code{\link{mLL}}, \code{\link{FULLmLL}} }
\examples{
   # Drosophila RNAi and Microarray Data from Boutros et al, 2002
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]

   # enumerate all possible models for 4 genes
   Sgenes = unique(colnames(D))
   models <- enumerate.models(Sgenes)  
   
   getRelevantEGenes(models[[64]], D, control=set.default.parameters(Sgenes, para=c(.13,.05), type="mLL"))

}

\keyword{models}
