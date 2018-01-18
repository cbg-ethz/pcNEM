\name{generateNetwork}
\alias{sampleRndNetwork}
\alias{sampleData}
\alias{sampleData.BN}

\title{Random networks and data sampling}
\description{1. Random network generation; 2. sampling of data from a given network topology}
\usage{
sampleRndNetwork(Sgenes, scaleFree=TRUE, gamma=2.5, maxOutDegree=length(Sgenes), maxInDegree=length(Sgenes), trans.close=TRUE, DAG=FALSE)

sampleData(Phi, m, prob=NULL, uninformative=0, type="binary", replicates=4, typeI.err=0.05, typeII.err=0.2, alpha=sample(seq(0.1,0.9,by=0.1),ncol(Phi),replace=TRUE), beta=sample(5:50,ncol(Phi),replace=TRUE), lambda=matrix(sample(seq(0.01,0.49,by=0.01),ncol(Phi)*2,replace=TRUE),ncol=2), meansH1=rep(0.5, ncol(Phi)), meansH0=rep(-0.5, ncol(Phi)), sdsH1=sample(seq(0.1,1,by=0.1),ncol(Phi),replace=TRUE), sdsH0=sample(seq(0.1,1,by=0.1),ncol(Phi),replace=TRUE))
}

\arguments{
\item{Sgenes}{character vector of S-genes}
\item{scaleFree}{should the network topology be scale free?}
\item{gamma}{for scale free networks: out-degrees of nodes are sampled from
  \eqn{\frac{1}{Z} * (0:maxOutDegree)^{-\gamma}}, where Z is a normalization
  factor}
 
\item{maxOutDegree}{maximal out-degree of nodes}
\item{maxInDegree}{maximal in-degree of nodes prior to transitive closure}
\item{trans.close}{Should the transitive closure of the graph be returned? Default: TRUE}
\item{DAG}{Should only DAGs be sampled? Default: FALSE}

\item{Phi}{adjacency matrix}
\item{m}{number of E-genes to sample}
\item{prob}{probability for each S-gene to get an E-gene attached}
\item{uninformative}{additional number of uninformative E-genes, i.e. E-genes carrying no information about the nested structure}
\item{type}{"binary" = binary data; "density" = log 'p-value' densities sampled from beta-uniform mixture model; "lodds" = log odds sampled from two normal distributions}
\item{replicates}{number of replicate measurements to simulate for binary data}
\item{typeI.err}{simulated type I error for binary data}
\item{typeII.err}{simulated type II error for binary data}
\item{alpha}{parameter for \eqn{Beta(\alpha,1)} distribution: one parameter per S-gene}
\item{beta}{parameter for \eqn{Beta(1,\beta)} distribution: one parameter per S-gene}
\item{lambda}{mixing coefficients for beta-uniform mixture model of the form:
  \eqn{\lambda_1 + \lambda_2*Beta(\alpha,1) + \lambda_3*Beta(1,\beta)}. There
  is a vector of 3 mixing coefficients per model and one model per S-gene.}
 
\item{meansH1}{normal distribution means of log odds ratios under the hypothesis of expecting an effect: one mean per S-gene}
\item{meansH0}{normal distribution means of log odds ratios under the null hypothesis: one mean per S-gene}
\item{sdsH1}{normal distribution standard deviations of log odds values under the hypothesis of expecting an effect: one sd per S-gene}
\item{sdsH0}{normal distribution standard deviations of log odds values under the null hypothesis: one sd per S-gene}
}
\details{
Random networks are generated as follows:
For each S-gene \eqn{S_{k}} we randomly choose the number \eqn{o} of outgoing
edges between 0 and maxOutDegree. This is either done uniform randomly or, if scale free networks are created, according to a power law distribution specified by gamma. We then select \eqn{o} S-genes having at most
maxInDegree ingoing edge and connected \eqn{S_{k}} to them. 

The function \code{sampleData} samples data from a given network topology as follows: We first attach E-genes to S-genes according to the probabilities prob (default: uniform). We then simulate knock-downs of the individual S-genes. For those E-genes, where no effects are expected, values are sampled from a null distribution, otherwise from an alternative distribution. In the simplest case we only sample binary data, where 1 indicates an effect an 0 no effect. Alternatively, we can sample log "p-value" densities according to a beta-uniform mixture model, where the null distribution is uniform and the alternative a mixture of two beta distributions. A third possibility is to sample log odds ratios, where alternative and null distribution are both normal.
}
\value{
For sampleRndNetwork an adjacency matrix, for sampleData a data matrix, for sampleData.BN a data matrix and a linking of effects to signals.
}

\author{Holger Froehlich, Cordula Zeller}


\seealso{\code{\link{getDensityMatrix}}}
\examples{
	Phi = sampleRndNetwork(paste("S",1:5,sep=""))
	D = sampleData(Phi, 100, type="density")$D	
  if(require(Rgraphviz)){
	  plot(as(transitive.reduction(Phi),"graphNEL"), main="original graph")
	  x11()
	  plot.nem(nem(D, control=set.default.parameters(unique(colnames(D)), type="CONTmLLBayes")), transitiveReduction=TRUE, SCC=FALSE, main    ="inferred graph")
  }
}
\keyword{graphs}
\keyword{models}
