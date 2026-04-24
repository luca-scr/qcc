\name{rules}
\alias{qccRules}
\alias{qccRulesViolatingWER1}
\alias{qccRulesViolatingWER2}
\alias{qccRulesViolatingWER3}
\alias{qccRulesViolatingWER4}
\alias{qccRulesViolatingNEL1}
\alias{qccRulesViolatingNEL2}
\alias{qccRulesViolatingNEL3}
\alias{qccRulesViolatingNEL4}
\alias{qccRulesViolatingNEL5}
\alias{qccRulesViolatingNEL6}
\alias{qccRulesViolatingNEL7}
\alias{qccRulesViolatingNEL8}

\title{Rules for Shewhart charts}

\description{
Functions which implement rules to signal out of control points in Shewhart charts.}

\usage{
qccRules(object, rules = object$rules, rule.set = object$rule.set)
}
\arguments{
  \item{object}{an object of class \code{'qcc'}.}
  \item{rules}{a vector of values specifying the rules used to declare a process out-of-control. Numeric values are interpreted within \code{rule.set}. The default rule set is \code{"western-electric"}, where possible values are:
  \describe{
  \item{1.}{One point plots outside 3-sigma control limits.}
  \item{2.}{Two of three consecutive points plot beyond a 2-sigma limit.}
  \item{3.}{Four of five consecutive points plot beyond a 1-sigma limit.}
  \item{4.}{Eight consecutive points plot on one side of the center line.}
  }
  These are often known as \emph{Western Electric rules}. With
  \code{rule.set = "nelson"}, values 1 through 8 request Nelson rules.
  }
  \item{rule.set}{a character string selecting the rule set used to interpret numeric \code{rules} values. Use \code{"western-electric"} for Western Electric rules or \code{"nelson"} for Nelson rules.}
}

\details{
The \code{qccRules()} function applies Western Electric rules by default.

The \code{qccRulesViolatingNEL1},
\code{qccRulesViolatingNEL2}, \dots, \code{qccRulesViolatingNEL8}
functions return the indices violating the corresponding Nelson rule:
\describe{
  \item{1.}{One point plots outside 3-sigma control limits.}
  \item{2.}{Nine points in a row plot on the same side of the center line.}
  \item{3.}{Six points in a row are steadily increasing or decreasing.}
  \item{4.}{Fourteen points in a row alternate up and down.}
  \item{5.}{Two of three consecutive points plot beyond a 2-sigma limit on the same side.}
  \item{6.}{Four of five consecutive points plot beyond a 1-sigma limit on the same side.}
  \item{7.}{Fifteen points in a row plot within 1 sigma of the center line.}
  \item{8.}{Eight points in a row plot outside 1 sigma on both sides of the center line.}
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
