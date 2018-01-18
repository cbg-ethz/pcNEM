\name{enumerate.models}
\alias{enumerate.models}

\title{Exhaustive enumeration of models}
\description{
  The function \code{enumerate.models} is used to create the model space for inference by exhaustive enumeration. 
  It computes a list of all transitively closed directed graphs on a given number of nodes.
}
\usage{
enumerate.models(x,name=NULL,trans.close=TRUE,verbose=TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{either the number of nodes or a vector of node names.}
  \item{name}{optionally the nodenames, if they are not provided in x}
  \item{trans.close}{should graphs be transitively closed?}
  \item{verbose}{if TRUE outputs number of (unique) models. Default: TRUE}
}
\details{
  The model space of Nested Effects Models consists of all transitively closed directed graphs. 
  The function \code{enumerate.models} creates them in three steps: (1.) build all directed graphs on \code{x} (or \code{length(x)}) nodes,
  (2.) transitively close each one of them, and (3.) remove redundant models to yield a unique set.
  So far, enumeration is limited to up to 5 nodes.
 
  I'm aware that this is inefficient! 
  It would be very desirable to enumerate the models directly (i.e. without creating all directed graphs as an intermediate step). 
}
\value{
a list of models. Each entry is a transitively closed adjacency matrix with unit main diagonal.
}
\author{Florian Markowetz}
\seealso{\code{\link{nem}}}
\examples{
enumerate.models(2)
enumerate.models(c("Anna","Bert"))
}
\keyword{graphs}
