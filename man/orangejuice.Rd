\name{orangejuice}
\alias{orangejuice}
\docType{data}
\title{Orange juice data}
\description{
Frozen orange juice concentrate is packed in 6-oz cardboard cans. These cans are formed on a machine by spinning them from cardboard stock and attaching a metal bottom panel. A can is then inspected to determine whether, when filled, the liquid could possible leak either on the side seam or around the bottom joint. If this occurs, a can is considered nonconforming.
The data were collected as 30 samples of 50 cans each at half-hour intervals over a three-shift period in which the machine was in continuous operation.
From sample 15 used a new batch of cardboard stock was punt into production.
Sample 23 was obtained when an inexperienced operator was temporarily assigned to the machine.
After the first 30 samples, a machine adjustment was made. Then further 24 samples were taken from the process.
}
\usage{data(orangejuice)}
\format{A data frame with 54 observations on the following 4 variables:
\describe{\item{sample}{sample id}
          \item{D}{number of defectives}
          \item{size}{sample sizes}
          \item{trial}{trial samples (TRUE/FALSE)}
         }
}
%\details{}
%\source{}
\references{
Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 152--155.
}
\examples{
data(orangejuice)
orangejuice$d <- orangejuice$D/orangejuice$size
attach(orangejuice)
summary(orangejuice)
boxplot(d ~ trial)
mark <- ifelse(trial, 1, 2)
plot(sample, d, type="b", col=mark, pch=mark)
detach(orangejuice)
}
\keyword{datasets}
