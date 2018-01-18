\name{plot.nem}
\alias{plot.nem}
\alias{plot.triples}
\alias{plot.pairwise}
\alias{plot.score}
\alias{plot.ModuleNetwork}
\alias{plot.nem.greedy}
\alias{plot.nem.greedyMAP}
\alias{plot.nem.bootstrap}
\alias{plot.nem.jackknife}
\alias{plot.nem.consensus}
\alias{plot.nem.BN}
\alias{plot.mc.eminem}
\alias{plot.dynoNEM}
\alias{plotnem}

\title{plot nested effect model}
\description{
  plot graph of nested effects model, the marginal likelihood distribution or the posterior position of the effected genes
}
\usage{
	\method{plot}{nem}(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, D=NULL, draw.lines=FALSE, palette="BlueRed",...)
}

\arguments{
 	\item{x}{nem object to plot}
	\item{what}{(i), "graph", (ii) "mLL" = likelihood distribution, (iii) "pos" = posterior position of effected genes}
	\item{remove.singletons}{remove unconnected nodes from the graph plot}
	\item{PDF}{output as PDF-file}
	\item{filename}{filename of PDF-file}
	\item{thresh}{if x has a real valued adjacency matrix (weight matrix),  don't plot edges with |weight| <= thresh}
	\item{transitiveReduction}{plot a transitively reduced graph} 
 	\item{plot.probs}{plot edge weights/probabilities. If regulation directions have been inferred (see \code{infer.edge.type}), upregulated edges are drawn red and downregulated edges blue. Edges, were no clear direction could be inferred, are drawn in black.} 
	\item{SCC}{plot the strongly connected components graph}
	\item{D}{Visualize the nested subset structure of the dataset via \code{plotEffects} along with the graph and show the linking of E-genes to S-genes in the dataset. Should only be used for small networks. Default: Just plot the graph}
	\item{draw.lines}{If the nested subset structure is shown, should additionally lines connecting S-genes and their associated E-genes be drawn? WARNING: For larger datasets than e.g. 5 S-genes this most probably does not work, because the nested subset structure picture then partially overlaps with the graph picture. Default: Do not draw these lines}
	\item{palette}{color palette to use: either 'BlueRed' (default) or 'Grey'}
	\item{...}{other arguments to be passed to the Rgraphviz plot function or to the graphics 'image' function.}
}
\value{
  none
}

\author{Florian Markowetz, Holger Froehlich}

\seealso{\code{\link{nem}}, \code{\link{plotEffects}}, \code{\link{infer.edge.type}}}

\keyword{models}
