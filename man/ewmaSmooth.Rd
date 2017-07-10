\name{ewmaSmooth}
\alias{ewmaSmooth}
\title{EWMA smoothing function}
\description{Compute Exponential Weighted Moving Average.}
\usage{
ewmaSmooth(x, y, lambda = 0.2, start, \dots)
}
\arguments{
\item{x}{a vector of x-values.}
\item{y}{a vector of y-values.}
\item{lambda}{the smoothing parameter.}
\item{start}{the starting value.}
\item{\dots}{additional arguments (currently not used).}
}
\details{
EWMA function smooths a series of data based on a moving average with weights which decay exponentially.

For each \eqn{y_t}{y_t} value the smoothed value is computed as
\deqn{z_t = \lambda y_t + (1-\lambda) z_{t-1}}
where \eqn{0 \le \lambda \le 1}{0 <= lambda <= 1} is the parameter which controls the weights applied.
}
\value{Returns a list with elements:
\item{x}{ordered x-values}
\item{y}{smoothed y-values}
\item{lambda}{the smoothing parameter}
\item{start}{the starting value}
}
\references{
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Wetherill, G.B. and Brown, D.W. (1991) \emph{Statistical Process Control}. New York: Chapman & Hall.
}
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
\seealso{\code{\link{qcc}}, \code{\link{cusum}}}
\examples{
x <- 1:50
y <- rnorm(50, sin(x/5), 0.5)
plot(x,y)
lines(ewmaSmooth(x,y,lambda=0.1), col="red")
}
\keyword{hplot}
