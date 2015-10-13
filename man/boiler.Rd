\name{boiler}
\alias{boiler}
\docType{data}
\title{Boiler temperature data}
\description{
Temperature readings from the eight configured burners on a boiler.
}
\usage{data(boiler)}
\format{
  A data frame with 25 observations on the following 8 variables:
  \describe{
    \item{t1}{temperature reading 1}
    \item{t2}{temperature reading 2}
    \item{t3}{temperature reading 3}
    \item{t4}{temperature reading 4}
    \item{t5}{temperature reading 5}
    \item{t6}{temperature reading 6}
    \item{t7}{temperature reading 7}
    \item{t8}{temperature reading 8}
  }
}
%\details{}
%\source{}
\references{Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM, p. 86.}
\examples{
data(boiler)
summary(boiler)
boxplot(boiler)
}
\keyword{datasets}
