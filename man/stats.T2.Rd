\name{stats.T2}
\alias{stats.T2}
\alias{limits.T2}

\title{Statistics used in computing and drawing the Hotelling T^2 chart for subgrouped data}
\description{These functions are used to compute statistics required by the \eqn{T^2} chart.}

\usage{
stats.T2(data, center = NULL, cov = NULL)

limits.T2(ngroups, size, nvars, conf) 
}

\arguments{
  \item{data}{the observed data values}
  \item{center}{a vector of values to use for center of input variables.}
  \item{cov}{a matrix of values to use for the covariance matrix of input variables.}
  \item{ngroups}{number of groups}
  \item{size}{sample size}
  \item{nvars}{number of variables}
  \item{conf}{confidence level (0 < \code{conf} < 1)}
}

%\details{}

\value{The function \code{stats.T2} returns a list with components: 
\item{statistics}{a vector of values for the \eqn{T^2} statistic}
\item{means}{a matrix of within group means for each variable}
\item{center}{sample/group center statistic}
\item{S}{covariance matrix}

The function \code{limits.T2} returns a list with components: 
\item{control}{control limits}
\item{prediction}{pred.limits}
}

\references{
Mason, R.L. and Young, J.C. (2002) \emph{Multivariate Statistical Process Control with Industrial Applications}, SIAM. \cr
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed. New York: John Wiley & Sons. \cr
Ryan, T. P. (2000), \emph{Statistical Methods for Quality Improvement}, 2nd ed. New York: John Wiley & Sons, Inc.
}

\author{Luca Scrucca}
%\note{}

\seealso{\code{\link{mqcc}}, \code{\link{stats.T2.single}}}
%\examples{}
\keyword{htest}
\keyword{hplot}
\keyword{multivariate}