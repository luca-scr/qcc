\name{describe}
\alias{describe}
\alias{print.describe}

\title{Descriptive statistics}

\description{This function compute and print simple descriptive statistics for a matrix or data frame.}

\usage{
describe <- function(data, \dots)
print.describe(x, digits = max(4, getOption("digits") - 3), \dots)
}

\arguments{
  \item{data}{a matrix or data.frame.}
  \item{\dots}{additional arguments to be passed to the generic function.}
  \item{digits}{the number of significant digits to use.}
}

\details{The function is able to deal with different types of variables, such as numeric, factor, and logical.}
%\value{}
%\references{}
\author{Luca Scrucca}
%\note{}
\seealso{\code{\link{summary}}}
\examples{
data(pistonrings)
describe(pistonrings)
by(pistonrings[,1:2], pistonrings$trial, describe)
}

