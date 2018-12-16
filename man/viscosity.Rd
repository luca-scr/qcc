\name{viscosity}
\alias{viscosity}
\docType{data}

\title{Viscosity for aircraft primer paint data}

\description{
"The viscosity of an aircraft primer paint is an important quality characteristic. The product is produced in batches, and because each batch takes several hours to produce, the production rate is too slow to allow for rational subgroups of size greater than one." (Montgomery, 2005, p. 232)
}

\usage{data(viscosity)}
\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{batch}{batch number}
    \item{viscosity}{viscosity measure}
    \item{trial}{preliminary sample indicator (TRUE/FALSE)}
  }
}
%\details{}
%\source{}
\references{
Montgomery, D.C. (2005) \emph{Introduction to Statistical Quality Control}, 5th ed, New York, John Wiley & Sons, pp. 232-235
}
\examples{
data(viscosity)
describe(viscosity, by = trial)
plot(viscosity ~ batch, data = viscosity, type = "o", 
     pch = ifelse(trial, 19, 1))
}
\keyword{datasets}
