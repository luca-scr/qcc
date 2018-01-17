\name{antifreeze}
\alias{antifreeze}

\title{Water content of antifreeze data}
\description{
Water content (in ppm) of batches of antifreeze.
}
\usage{data(antifreeze)}
\format{
  A data frame with observations on 34 successive batches of antifreeze:
  \describe{
    \item{water}{water content (in ppm)}
  }
}
%\details{}
%\source{}
\references{Wetherill, G.B. and Brown, D.W. (1991) Statistical Process Control. New York: Chapman & Hall, p. 120}
\examples{
data(antifreeze)
summary(antifreeze)
plot(antifreeze)
}
\keyword{datasets}
