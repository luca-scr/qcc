\name{qccGroups}
\alias{qccGroups}

\title{Grouping data based on a sample indicator}
\description{This function allows to easily group data to use as input to the \code{qcc()} function.}
\usage{
qccGroups(data, x, sample)
}
\arguments{
  \item{data}{a data frame (or a similar structure which can be coerced to be a \code{data.frame}) providing the observed data. If not provided, the input to the following arguments must be objects defined in the calling environment.}
  \item{x}{a name from \code{data} or a vector of observed data values.}
  \item{sample}{a name from \code{data} or a vector of sample indicators defining the rationale subgroups of data values.}
}

%\details{}
\value{
The function returns a matrix of suitable dimensions. If one or more group have fewer observations than others, \code{NA} values are used to fill empty values. 
}
%\references{}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{qcc}} }
\examples{
data(pistonrings)
# create a matrix of 40 samples made of 5 observations each
qccGroups(data = pistonrings, diameter, sample)
# remove some observations to get still a 40x5 matrix but filled with NAs 
qccGroups(data = pistonrings[-c(1,2,50,52,199),], diameter, sample)
}

