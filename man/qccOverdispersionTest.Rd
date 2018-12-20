\name{qccOverdispersionTest}
\alias{qccOverdispersionTest}

\title{Overdispersion test for binomial and poisson data}

\description{This function allows to test for overdispersed data in the binomial and poisson case.}

\usage{
qccOverdispersionTest(x, size, 
  type = ifelse(missing(size), "poisson", "binomial"))
}
\arguments{
  \item{x}{a vector of observed data values}
  \item{size}{for binomial data, a vector of sample sizes}
  \item{type}{a character string specifying the distribution for testing, either \code{"poisson"} or \code{"binomial"}. By default, if \code{size} is provided a binomial distribution is assumed, otherwise a poisson distribution.}
}
\details{
This very simple test amounts to compute the test statistic 
\deqn{D = s^2 / \sigma^2 \times (n - 1)}
where \eqn{s^2} is the observed variance, \eqn{\sigma^2} is the theoretical variance, and \eqn{n} is the number of observations. The test statistic is the compared to the critical value of a Chi-square distribution with \eqn{n-1} degrees of freedom.
}
\value{The function returns a matrix of results.}
\references{
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}, New York, Chapman and Hall, pp. 216--218
}
\author{Luca Scrucca}
%\note{}
%\seealso{}
\examples{
# data from Wetherill and Brown (1991) pp. 212--213, 216--218:
x <- c(12,11,18,11,10,16,9,11,14,15,11,9,10,13,12,
       8,12,13,10,12,13,16,12,18,16,10,16,10,12,14)
size <- rep(50, length(x))
qccOverdispersionTest(x, size)

x <- c(11,8,13,11,13,17,25,23,11,16,9,15,10,16,12,
       8,9,15,4,12,12,12,15,17,14,17,12,12,7,16)
qccOverdispersionTest(x)
}
\keyword{htest}

