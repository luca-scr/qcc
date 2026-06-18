

#' Water content of antifreeze data
#' 
#' Water content (in ppm) of batches of antifreeze.
#' 
#' 
#' @name antifreeze
#' @docType data
#' @format A data frame with observations on 34 successive batches of
#' antifreeze: \describe{ \item{water}{water content (in ppm)} }
#' @references Wetherill, G.B. and Brown, D.W. (1991) Statistical Process
#' Control. New York: Chapman & Hall, p. 120
#' @keywords datasets
#' @examples
#' 
#' data(antifreeze)
#' describe(antifreeze)
#' 
NULL





#' Boiler temperature data
#' 
#' Temperature readings from the eight configured burners on a boiler.
#' 
#' 
#' @name boiler
#' @docType data
#' @format A data frame with 25 observations on the following 8 variables:
#' \describe{ \item{t1}{temperature reading 1} \item{t2}{temperature reading 2}
#' \item{t3}{temperature reading 3} \item{t4}{temperature reading 4}
#' \item{t5}{temperature reading 5} \item{t6}{temperature reading 6}
#' \item{t7}{temperature reading 7} \item{t8}{temperature reading 8} }
#' @references Mason, R.L. and Young, J.C. (2002) \emph{Multivariate
#' Statistical Process Control with Industrial Applications}, SIAM, p. 86.
#' @keywords datasets
#' @examples
#' 
#' data(boiler)
#' describe(boiler)
#' boxplot(boiler)
#' 
NULL





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
#' 
#' @name circuit
#' @docType data
#' @format A data frame with 46 observations on the following 4 variables.
#' \describe{ \item{sample}{sample number} \item{x}{number of defectives in 100
#' printed circuit boards (inspection unit)} \item{size}{sample size}
#' \item{trial}{trial sample indicator (TRUE/FALSE)} }
#' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 173--175
#' @keywords datasets
#' @examples
#' 
#' data(circuit)
#' describe(circuit, by = trial)
#' boxplot(x/size ~ trial, data = circuit)
#' plot(x/size ~ sample, data = circuit, type="b")
#' 
NULL





#' Dyed cloth data
#' 
#' In a textile finishing plant, dyed cloth is inspected for the occurrence of
#' defects per 50 square meters. The data on ten rolls of cloth are presented.
#' 
#' 
#' @name dyedcloth
#' @docType data
#' @format A data frame with 10 observations on the following 2 variables.
#' \describe{ \item{x}{number of nonconformities per 50 square meters
#' (inspection units)} \item{size}{number of inspection units in roll (variable
#' sample size)} }
#' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 183--184
#' @keywords datasets
#' @examples
#' 
#' data(dyedcloth)
#' describe(dyedcloth)
#' dyedcloth  = transform(dyedcloth, sample = seq(nrow(dyedcloth)))
#' plot(x/size ~ sample, data = dyedcloth, type="b")
#' 
NULL





#' Orange juice data
#' 
#' Frozen orange juice concentrate is packed in 6-oz cardboard cans. These cans
#' are formed on a machine by spinning them from cardboard stock and attaching
#' a metal bottom panel. A can is then inspected to determine whether, when
#' filled, the liquid could possible leak either on the side seam or around the
#' bottom joint. If this occurs, a can is considered nonconforming. The data
#' were collected as 30 samples of 50 cans each at half-hour intervals over a
#' three-shift period in which the machine was in continuous operation. From
#' sample 15 used a new batch of cardboard stock was punt into production.
#' Sample 23 was obtained when an inexperienced operator was temporarily
#' assigned to the machine. After the first 30 samples, a machine adjustment
#' was made. Then further 24 samples were taken from the process.
#' 
#' 
#' @name orangejuice
#' @docType data
#' @format A data frame with 54 observations on the following 4 variables:
#' \describe{\item{sample}{sample id} \item{D}{number of defectives}
#' \item{size}{sample sizes} \item{trial}{trial samples (TRUE/FALSE)} }
#' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 152--155.
#' @keywords datasets
#' @examples
#' 
#' data(orangejuice)
#' orangejuice  = transform(orangejuice, d = D/size)
#' describe(orangejuice, by = trial)
#' boxplot(d ~ trial, data = orangejuice)
#' plot(d ~ sample, data = orangejuice, type = "b", pch = ifelse(trial, 1, 19))
#' 
NULL





#' Orange juice data -- Part 2
#' 
#' A full description of the problem is given in \code{\link{orangejuice}}. \cr
#' 
#' This dataset contains samples taken after the machine adjustment was made.
#' 
#' 
#' @name orangejuice2
#' @docType data
#' @format A data frame with 64 observations on the following 4 variables:
#' \describe{\item{sample}{sample id} \item{D}{number of defectives}
#' \item{size}{sample sizes} \item{trial}{trial samples (TRUE/FALSE)} }
#' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 155--159.
#' @keywords datasets
#' @examples
#' 
#' data(orangejuice2)
#' orangejuice2  = transform(orangejuice2, d = D/size)
#' describe(orangejuice2, by = trial)
#' boxplot(d ~ trial, data = orangejuice2)
#' plot(d ~ sample, data = orangejuice2, type = "b", pch = ifelse(trial, 1, 19))
#' 
NULL





#' Personal computer manufacturer data
#' 
#' A personal computer manufacturer counts the number of nonconformities per
#' unit on the final assembly line. He collects data on 20 samples of 5
#' computers each.
#' 
#' 
#' @name pcmanufact
#' @docType data
#' @format A data frame with 10 observations on the following 2 variables.
#' \describe{ \item{x}{number of nonconformities (inspection units)}
#' \item{size}{number of computers inspected} }
#' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 181--182
#' @keywords datasets
#' @examples
#' 
#' data(pcmanufact)
#' describe(pcmanufact)
#' pcmanufact  = transform(pcmanufact, sample = seq(nrow(pcmanufact)))
#' plot(x/size ~ sample, data = pcmanufact, type="b")
#' 
NULL





#' Piston rings data
#' 
#' Piston rings for an automotive engine are produced by a forging process. The
#' inside diameter of the rings manufactured by the process is measured on 25
#' samples, each of size 5, for the control phase I, when preliminary samples
#' from a process being considered 'in control' are used to construct control
#' charts. Then, further 15 samples, again each of size 5, are obtained for
#' phase II.
#' 
#' 
#' @name pistonrings
#' @docType data
#' @format A data frame with 200 observations on the following 3 variables.
#' \describe{ \item{diameter}{a numeric vector} \item{sample}{sample ID}
#' \item{trial}{preliminary sample indicator (TRUE/FALSE)} }
#' @references Montgomery, D.C. (1991) \emph{Introduction to Statistical
#' Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 206--213
#' @keywords datasets
#' @examples
#' 
#' data(pistonrings)
#' describe(pistonrings, by = trial)
#' boxplot(diameter ~ trial, data = pistonrings)
#' plot(diameter ~ sample, data = pistonrings, cex=0.7)
#' with(pistonrings, lines(tapply(diameter,sample,mean)))
#' 
NULL






#' Ryan's (2011) multivariate data for quality control
#' 
#' Multivariate data from Ryan (2011, Table 9.2).
#' 
#' 
#' @name RyanMultivar
#' @docType data
#' @format Multivariate data on 20 samples of size 4 for two variables.
#' \describe{ \item{RyanMultivar}{a list of two data frames, \code{X1} and
#' \code{X2}, one for each variable.} }
#' @references Ryan, T. P. (2011), \emph{Statistical Methods for Quality
#' Improvement}, 3rd ed. New York: John Wiley & Sons, Inc.
#' @keywords datasets
#' @examples
#' 
#' data(RyanMultivar)
#' str(RyanMultivar)
#' 
NULL





#' Viscosity for aircraft primer paint data
#' 
#' "The viscosity of an aircraft primer paint is an important quality
#' characteristic. The product is produced in batches, and because each batch
#' takes several hours to produce, the production rate is too slow to allow for
#' rational subgroups of size greater than one." (Montgomery, 2005, p. 232)
#' 
#' 
#' @name viscosity
#' @docType data
#' @format A data frame with 35 observations on the following 3 variables.
#' \describe{ \item{batch}{batch number} \item{viscosity}{viscosity measure}
#' \item{trial}{preliminary sample indicator (TRUE/FALSE)} }
#' @references Montgomery, D.C. (2005) \emph{Introduction to Statistical
#' Quality Control}, 5th ed, New York, John Wiley & Sons, pp. 232-235
#' @keywords datasets
#' @examples
#' 
#' data(viscosity)
#' describe(viscosity, by = trial)
#' plot(viscosity ~ batch, data = viscosity, type = "o", 
#'      pch = ifelse(trial, 19, 1))
#' 
NULL



