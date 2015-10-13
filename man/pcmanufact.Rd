\name{pcmanufact}
\alias{pcmanufact}
\docType{data}
\title{Personal computer manufacturer data}
\description{
A personal computer manufacturer counts the number of nonconformities per unit on the final assembly line. He collects data on 20 samples of 5 computers each.
}
\usage{data(pcmanufact)}
\format{
  A data frame with 10 observations on the following 2 variables.
  \describe{
    \item{x}{number of nonconformities (inspection units)}
    \item{size}{number of computers inspected}
  }
}
%\details{}
%\source{}
\references{Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 181--182}
\examples{
data(pcmanufact)
attach(pcmanufact)
summary(pcmanufact)
plot(x/size, type="b")
detach(pcmanufact)
}
\keyword{datasets}
