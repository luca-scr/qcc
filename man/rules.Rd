\name{rules}
\alias{qccRules}
\alias{qccRulesViolatingWER1}
\alias{qccRulesViolatingWER2}
\alias{qccRulesViolatingWER3}
\alias{qccRulesViolatingWER4}

\title{Rules for Shewhart charts}

\description{
Function which implements rules to signal out of control points in Shewhart charts.}

\usage{
qccRules(object, rules = object$rules)
}
\arguments{
  \item{object}{an object of class \code{'qcc'}.}
  \item{rules}{a vector of values specifying the rules used to declare a process out-of-control.
  Possible values are:
  \describe{
  \item{1.}{One point plots outside 3-sigma control limits.}
  \item{2.}{Two of three consecutive points plot beyond a 2-sigma limit.}
  \item{3.}{Four of five consecutive points plot beyond a 1-sigma limit.}
  \item{4.}{Eight consecutive points plot on one side of the center line.}
  }
  These are often known as \emph{Western Electric rules}.
  }
}

\value{
The function returns a vector of the same length as the statistics charted with value \code{NA} if no rule is violated, and numerical values corresponding to the violated rule. If more than one rule is violated, the value returned corresponds to the most serious violation.
}
%\references{ TODO: ~put references to the literature/web site here ~ }
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
%\examples{}
\keyword{hplot}
