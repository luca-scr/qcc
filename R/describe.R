#
# Descriptive statistics for a matrix or data frame.
#
# TODO: include also describe.by() ??

describe <- function(data, ...)
{
  if(is.matrix(data) | is.data.frame(data))
    { varnames <- dimnames(data)[[2]] }
  else
    { varnames <- deparse(substitute(data))
      data <- as.matrix(data) }
  nvar <- ncol(data)
  obj <- vector(mode="list", length=nvar)
  names(obj) <- varnames
  type <- rep(NA, nvar)

  opt.warn <- options("warn")  # save default warning option
  options(warn=-1)             # and suppress warnings
  for(i in seq(nvar))
     { 
       x <- data[,i]
       if(is.factor(x) | typeof(x) == "character" | typeof(x) == "logical")
         { type[i] <- "factor"
           out <- summary.factor(x)
           obj[[i]] <- out  
       }
       else if(any(class(x) == "POSIXt"))
         { type[i] <- "numeric"
           out <- summary(x)
           obj[[i]] <- out  
       }
       else
         { type[i] <- "numeric"
           n.miss <- sum(is.na(x))
           x <- na.omit(x)
           out <- c(length(x), mean(x), sd(x), fivenum(x))
           names(out) <- c("Obs", "Mean", "Std.Dev.", "Min", "Q1", "Median", "Q3", "Max")
           obj[[i]] <- out
         }
  }
  obj <- list(describe = obj, type = type)
  class(obj) <- "describe"
  options(warn=opt.warn$warn)
  return(obj)
}

print.describe <- function(x, digits = max(4, getOption("digits") - 3), ...)
{
  descr <- x$describe
  isNum <- (x$type == "numeric")
  
  if(sum(isNum) > 0)
  { 
    out1 <- do.call("rbind",descr[isNum])
    print(out1, digits = digits)
  }

  if(sum(!isNum) > 0)
  {
    out2 <- descr[!isNum]
    for(j in seq(out2))
    { 
      tab <- out2[[j]]
      tab <- cbind("Freq." = tab, "Rel.Freq." = tab/sum(tab), 
                   "Cum.Rel.Freq." = cumsum(tab)/sum(tab))
      names(dimnames(tab)) <- c(names(out2)[j], "")
      print(tab, digits = digits)
    }
  }
  
  invisible()
}
