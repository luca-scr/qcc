#-----------------------------------------------------------------------------#
#                                                                             # 
#  Example showing how to extend the package by defining a new control chart  #
#                                                                             #
# Reference:                                                                  #
# Scrucca, L. (2004). qcc: an R package for quality control charting and      #
#   statistical process control. R News, 4/1, 11-17.                          #
#-----------------------------------------------------------------------------#

# Standardized p chart defined as type = "p.std"

# Compute group statistics and center
stats.p.std <- function(data, sizes)
{
  data <- as.vector(data)
  sizes <- as.vector(sizes)
  pbar <- sum(data)/sum(sizes)
  z <- (data/sizes - pbar)/sqrt(pbar*(1-pbar)/sizes)
  list(statistics = z, center = 0)
}

# Compute within-group standard deviation
sd.p.std <- function(data, sizes, ...) { return(1) }

# Compute control limits
limits.p.std <- function(center, std.dev, sizes, conf)
{
  if(conf >= 1) 
    { lcl <- -conf
      ucl <- +conf 
  }
  else
    { if(conf > 0 & conf < 1)
        { nsigmas <- qnorm(1 - (1 - conf)/2)
          lcl <- -nsigmas
          ucl <- +nsigmas }
      else stop("invalid 'conf' argument.") 
  }
  limits <- matrix(c(lcl, ucl), ncol = 2)
  rownames(limits) <- rep("", length = nrow(limits))
  colnames(limits) <- c("LCL", "UCL")
  return(limits)
}

# Example with simulated data

# set unequal sample sizes
n <- c(rep(50,5), rep(100,5), rep(25, 5))
# generate randomly the number of successes
x <- rbinom(length(n), n, 0.2)
# plot the control chart with variable limits
summary(qcc(x, type="p", size=n))
# plot the standardized control chart
summary(qcc(x, type="p.std", size=n))

