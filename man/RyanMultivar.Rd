\name{RyanMultivar}
\alias{RyanMultivar}
\docType{data}

\title{Ryan's (2011) multivariate data for quality control}

\description{
Multivariate data from Ryan (2011, Table 9.2).
}

\usage{data(RyanMultivar)}
\format{
  Multivariate data on 20 samples of size 4 for two variables.
  \describe{
    \item{RyanMultivar}{a list of two data frames, \code{X1} and \code{X2}, 
          one for each variable.}
  }
}
%\details{}
%\source{}
\references{
Ryan, T. P. (2011), \emph{Statistical Methods for Quality Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.
}
\examples{
data(RyanMultivar)
str(RyanMultivar)
}
\keyword{datasets}
