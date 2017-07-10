\name{shewhart.rules}
\alias{shewhart.rules}
\alias{beyond.limits}
\alias{violating.runs}
\title{Functions specifying rules for Shewhart charts}
\description{
These functions are used to signal out of control points in Shewhart charts.}
\usage{
shewhart.rules(object, limits = object$limits, run.length = qcc.options("run.length"))
beyond.limits(object, limits = object$limits)
violating.runs(object, run.length = qcc.options("run.length"))
}
\arguments{
  \item{object}{an object of class 'qcc'.}
  \item{limits}{control limits}
  \item{run.length}{the maximum value of a run before to signal a
  point as out of control.}
}
\details{
The function \code{shewhart.rules} simply calls the \code{beyond.limits} and \code{violating.runs} functions which actually do the real calculations.
}
\value{
The \code{shewhart.rules} function returns a list with components:
\item{beyond.limits}{the indices of points beyond control limits.}
\item{violating.runs}{the indices of points violating runs.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Luca Scrucca}
%\note{ ~~further notes~~ }
%\examples{}
\keyword{hplot}
