\name{qcc.groups}
\alias{qcc.groups}

\title{Grouping data based on a sample indicator}
\description{This function allows to easily group data to use as input to the 'qcc' function.}
\usage{
qcc.groups(data, sample)
}
\arguments{
  \item{data}{the observed data values}
  \item{sample}{the sample indicators for the data values}
}
%\details{}
\value{
The function returns a matrix of suitable dimensions. If one or more group have few observations than others, \code{NA} values are appended.
}
%\references{}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
\examples{
data(pistonrings)
attach(pistonrings)
# 40 sample of 5 obs each
qcc.groups(diameter, sample)
# some obs are removed, the result is still a 40x5 matrix but with NAs added
qcc.groups(diameter[-c(1,2,50,52, 199)], sample[-c(1,2,50,52, 199)])
}
\keyword{manip}

