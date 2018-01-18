\name{plotEffects}
\alias{plotEffects}

\title{Plots data according to a phenotypic hierarchy}
\description{
  \code{plotEffects} visualizes the subset structure in the data by reordering rows and columns according to the topological order given by a phenotypic hierarchy.
}
\usage{
plotEffects(D,nem,border=TRUE,legend=TRUE,order=NULL,orderSCC=TRUE,palette="BlueRed",...)
}

\arguments{
  \item{D}{data matrix}
  \item{nem}{phenotypic hierarchy (object of class 'score' or 'pairwise')}
  \item{border}{draw red lines to indicate gene-specific effect reporters. Default: TRUE}
  \item{legend}{plot a legend. Default: TRUE}
  \item{order}{pre-define an order of the S-genes instead of the topological order to visualize the subset structure. Default: Use topological order.}
  \item{orderSCC}{Is the pre-defined order given on strongly connected components rather than on individual nodes?}
  \item{palette}{color palette to use: either 'BlueRed' (default) or 'Grey'}
  \item{...}{additional parameters for the graphics function 'image'}
}
\details{
  The experiments in the columns are reordered according to the topological order given by a phenotypic hierarchy.
  The effect reporters in the rows are grouped together by their position in the hierarchy.
  The groups are then arranged by topological order.
  Within each group the rows are hierarchically clustered.
  }
\value{
	ordering of the E-genes according to the hierarchy (vector of indices)
}
\note{
	This function was formerly named \code{plot.effects}. This naming is not possible any more, since S3 classes were used for the function \code{plot.nem}.
}

\author{Florian Markowetz, Holger Froehlich}

\examples{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   res <- nem(D,control=set.default.parameters(unique(colnames(D)), para=c(.13,.05)))
   if(require(Rgraphviz))
    plotEffects(D,res)
}
\keyword{graphs}
