\name{dyedcloth}
\alias{dyedcloth}
\docType{data}
\title{Dyed cloth data}
\description{
In a textile finishing plant, dyed cloth is inspected for the occurrence of defects per 50 square meters. The data on ten rolls of cloth are presented.
}
\usage{data(dyedcloth)}
\format{
  A data frame with 10 observations on the following 2 variables.
  \describe{
    \item{x}{number of nonconformities per 50 square meters (inspection units)}
    \item{size}{number of inspection units in roll (variable sample size)}
  }
}
%\details{}
%\source{}
\references{Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 183--184}
\examples{
data(dyedcloth)
attach(dyedcloth)
summary(dyedcloth)
plot(x/size, type="b")
detach(dyedcloth)
}
\keyword{datasets}
