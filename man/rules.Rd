\name{rules}
\alias{qccRules}
\alias{shewhart.rules}
\alias{beyond.limits}
\alias{violating.runs}
\alias{shewhartRules}
\alias{qccRulesViolatingWER1}
\alias{qccRulesViolatingWER2}
\alias{qccRulesViolatingWER3}
\alias{qccRulesViolatingWER4}

\title{Functions specifying rules for Shewhart charts}
\description{
These functions are used to signal out of control points in Shewhart charts.}

\usage{

shewhart.rules(object, limits = object$limits, run.length = qcc.options("run.length"))

beyond.limits(object, limits = object$limits)

violating.runs(object, run.length = qcc.options("run.length"))

qccRules(object, rules = object$rules)
}
\arguments{
  \item{object}{an object of class \code{'qcc'}.}
  \item{limits}{control limits}
  \item{run.length}{the maximum value of a run before to signal a
  point as out of control.}
  \item{rules}{a vector of values specifying the rules used to declare a process out-of-control.
  Possible values are:
  \describe{
  \item{1.}{One point plots outside 3-sigma control limits.}
  \item{2.}{Two of three consecutive points plot beyond a 2-sigma limit.}
  \item{3.}{Four of five consecutive points plot beyond a 1-sigma limit.}
  \item{4.}{Eight consecutive points plot on one side of the center line.}
  }
  These are often known as Western Electric rules.
  }
}
\details{
TODO: to add recently added function for WER and qccRules. 
TODO: To review all
The function \code{shewhart.rules} simply calls the \code{beyond.limits} and \code{violating.runs} functions which actually do the real calculations.
}
\value{
The \code{shewhart.rules} function returns a list with components:
\item{beyond.limits}{the indices of points beyond control limits.}
\item{violating.runs}{the indices of points violating runs.}
}
%\references{ TODO: ~put references to the literature/web site here ~ }
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
%\examples{}
\keyword{hplot}
