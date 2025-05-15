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
orangejuice2  = transform(orangejuice2, d = D/size)
describe(orangejuice2, by = trial)
boxplot(d ~ trial, data = orangejuice2)
plot(d ~ sample, data = orangejuice2, type = "b", pch = ifelse(trial, 1, 19))
}
\keyword{datasets}
