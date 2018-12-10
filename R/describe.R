#
# Descriptive statistics for a matrix or data frame
#

describe <- function(data, by, detail = FALSE, ...)
{
  data_name <- deparse(substitute(data))
  if(!is.data.frame(data)) 
    data <- as.data.frame(data)
  varnames <- colnames(data)
    
  if(!missing(by))
  {
    # TODO: check if ok when data is a matrix/data.frame
    by_name <- deparse(substitute(by))
    if(is.null(data[[by_name]]))
      stop(cat(by_name, "not available in data.frame", data_name))
    by <- as.factor(data[[by_name]])
    x <- split(data[setdiff(varnames, by_name)], by)
    out <- vector("list", length = nlevels(by))
    for(i in seq(nlevels(by)))
    {
      out[[i]] <- describe(x[[i]], detail = detail)
    }
    names(out) <- levels(by)
    out$by <- by_name
    class(out) <- c("describe")
    return(out)
  }
  
  nvar <- length(varnames)
  obj <- vector(mode="list", length=nvar)
  names(obj) <- if(nvar > 1) varnames else data_name
  type <- rep(NA, nvar)

  opt.warn <- options("warn")  # save default warning option
  options(warn=-1)             # and suppress warnings
  #
  skewness <- function(x) mean((x - mean(x))^3)/sd(x)^3
  kurtosis <- function(x) mean((x - mean(x))^4)/sd(x)^4 - 3
  #
  for(j in seq(nvar))
  { 
    x <- data[,j]
    if(is.factor(x) | typeof(x) == "character" | typeof(x) == "logical")
    { 
      type[j] <- "factor"
      if(detail)
      {
        out <- summary(as.factor(x))
        out <- cbind("Freq." = out, "Percent" = out/sum(out)*100, 
                     "Cum.Percent" = cumsum(out/sum(out)*100))
      } else
      {
        out <- summary(as.factor(x))
      }
      obj[[j]] <- out  
    }
    else if(any(class(x) == "POSIXt"))
    { 
      type[j] <- "numeric"
      out <- summary(x)
      obj[[j]] <- out  
    }
    else
    { 
      type[j] <- "numeric"
      n.miss <- sum(is.na(x))
      x <- na.omit(x)
      if(detail)
      {
        out <- c(length(x), n.miss, mean(x), sd(x), fivenum(x),
                 skewness(x), kurtosis(x))
        names(out) <- c("Obs", "NAs", "Mean", "Std.Dev.", 
                        "Min", "Q1", "Median", "Q3", "Max",
                        "Skewness", "Kurtosis")
      } else
      {
        out <- c(length(x), mean(x), sd(x), min(x), max(x))
        names(out) <- c("Obs", "Mean", "Std.Dev.", "Min", "Max")
      }
      obj[[j]] <- out
    }
  }
  obj <- list(name = data_name, describe = obj, type = type)
  class(obj) <- "describe"
  options(warn=opt.warn$warn)
  return(obj)
}

print.describe <- function(x, digits = getOption("digits") - 3, ...)
{

  if(!is.null(x$by))
  {
    by <- which(sapply(x, class) == "describe")
    for(i in by)
    {
      cat(cli::rule(left = paste(x$by, "=", names(x)[i])), "\n\n")
      print(x[[i]])
      if(i < length(by)) cat("\n")
    }
    return(invisible())  
  }
  
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
      if(is.vector(out2[[j]]))
      {
        out2[[j]] <- do.call("rbind", out2[j])
        cat("\n")
        print(out2[[j]], digits = digits)
      } else
      {
        names(dimnames(out2[[j]])) <- c(names(out2)[j], "")
        print(out2[[j]], digits = digits)
      }
    }
  }
  
  invisible()
}
