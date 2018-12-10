\name{describe}
\alias{describe}
\alias{print.describe}

\title{Descriptive statistics}

\description{Compute and print descriptive statistics for a matrix or data frame, eventually conditioning on another variable.}

\usage{
describe(data, by, detail = FALSE, \dots)
\method{print}{describe}(x, digits = getOption("digits") - 3, \dots)
}

\arguments{
  \item{data}{a matrix or data.frame.}
  \item{by}{a factor or character vector specifying the conditioning variable.}
  \item{detail}{a logical specifying if detailed statistics should be provided.}
  \item{\dots}{additional arguments to be passed to the generic function.}
  \item{x}{an object of class 'describe'.}
  \item{digits}{the number of significant digits to use.}
}

\details{This function provides simple descriptive statistics for different types of variables, such as numeric, factor, and logical, conditioning or not on another variable.}

%\value{}
%\references{}
\author{Luca Scrucca}
%\note{}

\seealso{\code{\link{summary}}, \code{\link{by}}}
\examples{
data(warpbreaks)
describe(warpbreaks)
describe(warpbreaks, detail = TRUE)
describe(warpbreaks, by = wool)
}

