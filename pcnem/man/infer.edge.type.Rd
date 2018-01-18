\name{infer.edge.type}
\alias{infer.edge.type}

\title{Infer regulation direction for each edge}
\description{
	The method infers edge types (up-regulation, down-regulation) for a given nem model. 
	Direct approach: For an edge a->b the method looks, whether b is up- or down-regulated in a knock-down of a.
	
	Indirect approach: For an edge a->b the method looks at the fraction of E-genes attached to b (including b itself), which are up- or down-regulated in a knock-down of a. If significantly more genes are down-regulated than up-regulated, the edge a->b is assumed to be an activation. Likewise, if significantly more genes are up-regulated than down-regulated, a->b is assumed to be an inhibition. If there is no significant difference in up- and down-regulated edges, a->b does not have a specified type.
}
\usage{
infer.edge.type(x, logFC, alpha=0.05, adj.method="BY", method=c("direct", "indirect"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{  
  \item{x}{nem object}    
  \item{logFC}{matrix with fold changes. The rownames of this matrix should correspond to the rownames of the data matrix, which was used to infer the nem model.}
  \item{alpha}{p-value cutoff}
  \item{adj.method}{multiple testing correction method. Default: Benjamini-Yekutieli}  
  \item{method}{default: "direct"}
}
\details{
	Significance in case of the indirect method is calculated using a two-tailed binomial test with null hypothesis p=0.5.
}
\value{
  Modified nem object. Each edge in the nem graph now has a "weight" and a "label" attribute. The label attribute corresponds to the original value in the adjacency matrix. The weight attribute encodes up- and down-regulation in the following way: value 2 means up-regulation, value -1 down-regulation and value 1 an unknown effect.
}
\author{Holger Froehlich}
\seealso{\code{\link{binom.test}}}

\examples{
   	data("BoutrosRNAi2002") 
	D <- BoutrosRNAiDiscrete[,9:16]   	
	result = nem(D, control=set.default.parameters(unique(colnames(D)), para=c(0.13,0.05)))
   	resEdgeInf = infer.edge.type(result, BoutrosRNAiLogFC)
  if(require(Rgraphviz))
   	plot.nem(resEdgeInf)
}

\keyword{models}
