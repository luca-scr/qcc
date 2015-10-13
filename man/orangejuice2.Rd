\name{orangejuice2}
\alias{orangejuice2}
\docType{data}
\title{Orange juice data --  Part 2}
\description{
A full description of the problem is given in \code{\link{orangejuice}}. \cr

This dataset contains samples taken after the machine adjustment was made. 
}
\usage{data(orangejuice)}
\format{
A data frame with 64 observations on the following 4 variables:
\describe{\item{sample}{sample id}
          \item{D}{number of defectives}
          \item{size}{sample sizes}
          \item{trial}{trial samples (TRUE/FALSE)}
         }
}
%\details{}
%\source{}
\references{
Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 155--159.
}
\examples{
data(orangejuice2)
orangejuice2$d <- orangejuice2$D/orangejuice2$size
attach(orangejuice2)
summary(orangejuice2)
boxplot(d ~ trial)
mark <- ifelse(trial, 1, 2)
plot(sample, d, type="b", col=mark, pch=mark)
detach(orangejuice2)
}
\keyword{datasets}
