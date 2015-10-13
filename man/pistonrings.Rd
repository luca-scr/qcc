\name{pistonrings}
\alias{pistonrings}
\docType{data}
\title{Piston rings data}
\description{
Piston rings for an automotive engine are produced by a forging process. The inside diameter of the rings manufactured by the process is measured on 25 samples, each of size 5, for the control phase I, when preliminary samples from a process being considered 'in control' are used to construct control charts. Then, further 15 samples, again each of size 5, are obtained for phase II.
}
\usage{data(pistonrings)}
\format{
  A data frame with 200 observations on the following 3 variables.
  \describe{
    \item{diameter}{a numeric vector}
    \item{sample}{sample ID}
    \item{trial}{preliminary sample indicator (TRUE/FALSE)}
  }
}
%\details{}
%\source{}
\references{Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 206--213}
\examples{
data(pistonrings)
attach(pistonrings)
summary(pistonrings)
boxplot(diameter ~ sample)
plot(sample, diameter, cex=0.7)
lines(tapply(diameter,sample,mean))
detach(pistonrings)
}
\keyword{datasets}
