[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/qcc)](https://cran.r-project.org/package=qcc)
[![](http://cranlogs.r-pkg.org/badges/qcc)](https://cran.r-project.org/package=qcc)
[![](http://cranlogs.r-pkg.org/badges/grand-total/qcc)](https://cran.r-project.org/package=qcc)

# qcc: An R package for quality control charting and statistical process control

This package provides quality control tools for statistical process control. 

* Shewhart quality control charts for continuous, attribute and count data. 
* Cusum and EWMA charts. 
* Operating characteristic curves. 
* Process capability analysis. 
* Pareto chart and cause-and-effect chart. 
* Multivariate control charts.

#### Author/Maintainer

* Luca Scrucca

#### Contributors

* Greg Snow
* Peter Bloomfield

## Installation

To install the package type the following:

```
install.packages("qcc")
library(qcc)
```

Or you can install the development version from GitHub:

```
library(devtools)
install_github("luca-scr/qcc")
library(qcc)
```

## How to Use This Package

See the paper 

> Scrucca, L. (2004) qcc: an R package for quality control charting and statistical process control. *R News* 4/1, 11-17.

and the vignette 

> A quick tour of qcc

Both documents are available in the **Vignettes and other documentation** section of the help.