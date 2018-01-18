\name{BFSlevel}
\alias{BFSlevel}

\title{Build (generalized) hierarchy by Breath-First Search}
\description{
  \code{BFSlevel} builds a (generalized) hierarchy by Breath-First Search as described in (Yu and Gerstein, 2006)
}
\usage{
BFSlevel(g,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{graphNEL object}
  \item{verbose}{Default: TRUE}
}

\details{
Haiyuan Yu and Mark Gerstein: Genomic analysis of the hierarchical structure of regulatory networks, PNAS 103(40):14724-14731, 2006 
}

\value{
  \item{level}{vector of levels for each node}

}

\author{Florian Markowetz}


%\seealso{}
\examples{
 ## bla
}
\keyword{graphs}
