\name{qccGroups}
\alias{qccGroups}

\title{Grouping data based on a sample indicator}
\description{This function allows to easily group data to use as input to the \code{qcc()} function.}
\usage{
qccGroups(x, sample, data)
}
\arguments{
  \item{x}{a vector of observed data values}
  \item{sample}{a vector of sample indicators for the data values}
  \item{data}{a data frame providing the observed data values \code{x} and the indicators \code{sample}. If not provided, the previous arguments must be objects defined in the calling environment.}
}
%\details{}
\value{
The function returns a matrix of suitable dimensions. If one or more group have fewer observations than others, \code{NA} values are appended.
}
%\references{}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
\examples{
data(pistonrings)
# create a matrix of 40 samples made of 5 observations each
qccGroups(diameter, sample, data = pistonrings)
# if some observations are removed, the result is still a 40x5 matrix but 
# with NAs added
qccGroups(diameter, sample, data = pistonrings[-c(1,2,50,52,199),])
}

