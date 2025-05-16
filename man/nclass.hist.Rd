\name{nclass.hist}
\alias{nclass.hist}

\title{Histogram number of classes/bins}

\usage{
nclass.hist(x, ...)
}

\arguments{
\item{x}{a vector of data values.}

\item{...}{additional arguments to be passed to low level functions.}
}

\description{
Computes the optimal number of classes/bins for an histogram as the
maximum between the Sturges and Freedman-Diaconis (FD) estimators.
For small datasets the Sturges value will usually be chosen, while
larger datasets will usually default to FD. Avoids the overly
conservative behaviour of FD and Sturges for small and large
datasets, respectively.
This is the default option in \code{numpy.histogram_bin_edges} available
in Python.
}

\value{The value of suggested number of classes/bins.}

\author{Luca Scrucca}

\seealso{\code{\link[grDevices]{nclass.FD}}, 
         \code{\link[grDevices]{nclass.Sturges}}}

\examples{
set.seed(1)
x <- stats::rnorm(111)
nclass.hist(x)
x <- stats::rnorm(1111)
nclass.hist(x)
}
