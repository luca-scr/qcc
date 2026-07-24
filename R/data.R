# TODO: data: data.frame -> numeric vector

#' Water content of antifreeze data
#'
#' Water content in 34 successive batches of antifreeze.
#'
#' @format A data frame with 34 rows and 1 column:
#' \describe{
#'   \item{water}{water content, in ppm}
#' }
#'
#' @source Wetherill, G.B. and Brown, D.W. (1991) Statistical Process
#' Control. New York: Chapman & Hall, p. 120
#'
#' @family datasets in `qcc` package
"antifreeze"




# TODO: docs: add units of measurement to each column
# TODO: docs: elaborate more in description

#' Boiler temperature data
#'
#' Temperature readings from the eight configured burners on a boiler.
#'
#' @format A data frame with 25 observations on the following 8 variables:
#' \describe{
#'   \item{t1}{temperature reading 1}
#'   \item{t2}{temperature reading 2}
#'   \item{t3}{temperature reading 3}
#'   \item{t4}{temperature reading 4}
#'   \item{t5}{temperature reading 5}
#'   \item{t6}{temperature reading 6}
#'   \item{t7}{temperature reading 7}
#'   \item{t8}{temperature reading 8}
#' }
#'
#' @source Mason, R.L. and Young, J.C. (2002) \emph{Multivariate
#' Statistical Process Control with Industrial Applications}, SIAM, p. 86.
#'
#' @family datasets in `qcc` package
"boiler"





#' Circuit boards data
#'
#' Number of nonconformities observed in 26 successive samples of 100 printed
#' circuit boards. Sample 6 and 20 are outside the control limits. Sample 6 was
#' examined by a new inspector and he did not recognize several type of
#' nonconformities that could have been present. Furthermore, the unusually
#' large number of nonconformities in sample 20 resulted from a temperature
#' control problem in the wave soldering machine, which was subsequently
#' repaired. The last 20 samples are further samples collected on inspection
#' units (each formed by 100 boards).
#'
#' @format A data frame with 46 observations on the following 4 variables:
#' \describe{
#'   \item{sample}{sample number}
#'   \item{x}{number of defectives in 100 printed circuit boards (inspection unit)}
#'   \item{size}{sample size}
#'   \item{trial}{trial sample indicator (TRUE/FALSE)}
#' }
#'
#' @source Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 173--175
#'
#' @family datasets in `qcc` package
"circuit"





#' Dyed cloth data
#'
#' In a textile finishing plant, dyed cloth is inspected for the occurrence of
#' defects per 50 square meters. The data on ten rolls of cloth are presented.
#'
#' @format A data frame with 10 observations on the following 2 variables:
#' \describe{
#'   \item{x}{number of nonconformities per 50 square meters (inspection units)}
#'   \item{size}{number of inspection units in roll (variable sample size)}
#' }
#'
#' @source Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 183--184
#'
#' @family datasets in `qcc` package
"dyedcloth"





#' Orange juice data
#'
#' Frozen orange juice concentrate is packaged in 6-ounce cardboard cans. A
#' machine forms each can by spinning it from cardboard stock and attaching a
#' metal bottom panel. A can is classified as nonconforming if an inspection
#' indicates that it could leak when filled, either along the side seam or
#' around the bottom joint.
#'
#' Thirty samples of 50 cans each were collected at half-hour intervals while
#' the machine operated continuously over three shifts. A new batch of cardboard
#' stock was introduced beginning with sample 15, and an inexperienced operator
#' was temporarily assigned to the machine for sample 23. After these 30
#' samples, the machine was adjusted and 24 additional samples were collected.
#' The `orangejuice2` dataset contains samples collected after this adjustment.
#'
#' @format `orangejuice` is a data frame with 54 observations and `orangejuice2`
#' is a data frame with 64 observations. Both contain the following 4 variables:
#' \describe{
#'   \item{sample}{sample id}
#'   \item{D}{number of defectives}
#'   \item{size}{sample sizes}
#'   \item{trial}{trial samples (TRUE/FALSE)}
#' }
#'
#' @source Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 152--159.
#'
#' @family datasets in `qcc` package
"orangejuice"

#' @rdname orangejuice
#' @format NULL
"orangejuice2"





#' Personal computer manufacturer data
#'
#' A personal computer manufacturer counts the number of nonconformities per
#' unit on the final assembly line. He collects data on 20 samples of 5
#' computers each.
#'
#' @format A data frame with 10 observations on the following 2 variables:
#' \describe{
#'   \item{x}{number of nonconformities (inspection units)}
#'   \item{size}{number of computers inspected}
#' }
#'
#' @source Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 181--182
#'
#' @family datasets in `qcc` package
"pcmanufact"




# TODO: docs: add units of measurement to diameter

#' Piston rings data
#'
#' Piston rings for an automotive engine are produced by a forging process. The
#' inside diameter of the rings manufactured by the process is measured on 25
#' samples, each of size 5, for the control phase I, when preliminary samples
#' from a process being considered 'in control' are used to construct control
#' charts. Then, further 15 samples, again each of size 5, are obtained for
#' phase II.
#'
#' @format A data frame with 200 observations on the following 3 variables:
#' \describe{
#'   \item{diameter}{a numeric vector}
#'   \item{sample}{sample ID}
#'   \item{trial}{preliminary sample indicator (TRUE/FALSE)}
#' }
#'
#' @source Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 206--213
#'
#' @family datasets in `qcc` package
"pistonrings"




# TODO: docs: rewrite the format section

#' Ryan's (2011) multivariate data for quality control
#'
#' Multivariate data from Ryan (2011, Table 9.2).
#'
#' @format Multivariate data on 20 samples of size 4 for two variables:
#' \describe{
#' \item{RyanMultivar}{a list of two data frames, \code{X1} and
#' \code{X2}, one for each variable.}
#' }
#'
#' @source Ryan, T. P. (2011), \emph{Statistical Methods for Quality
#' Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.
#'
#' @family datasets in `qcc` package
"RyanMultivar"




# TODO: docs: add units of measurement to viscosity

#' Viscosity for aircraft primer paint data
#'
#' "The viscosity of an aircraft primer paint is an important quality
#' characteristic. The product is produced in batches, and because each batch
#' takes several hours to produce, the production rate is too slow to allow for
#' rational subgroups of size greater than one." (Montgomery, 2005, p. 232)
#'
#' @format A data frame with 35 observations on the following 3 variables.
#' \describe{
#'   \item{batch}{batch number}
#'   \item{viscosity}{viscosity measure}
#'   \item{trial}{preliminary sample indicator (TRUE/FALSE)}
#' }
#'
#' @source Montgomery, D.C. (2005) \emph{Introduction to Statistical
#' Quality Control}, 5th ed, New York, John Wiley & Sons, pp. 232-235
#'
#' @family datasets in `qcc` package
"viscosity"
