#-------------------------------------------------------------------#
#                                                                   #
#    Process Capability Analysis                                    #
#                                                                   #
#-------------------------------------------------------------------#

processCapability <- function(object, spec.limits, target, 
                               std.dev, nsigmas, confidence.level = 0.95,
                               plot = TRUE, ...)
{
# Computes process capability indices for a qcc object of type "xbar" 
# and plot the histogram

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")
  if (!(object$type=="xbar" | object$type=="xbar.one"))
     stop("Process Capability Analysis only available for charts type \"xbar\" and \"xbar.one\" charts")

  x <- as.vector(object$data)
  x <- x[!is.na(x)]
  center <- object$center
  if(missing(std.dev))
    std.dev <- object$std.dev
  n <- length(x)

  if (missing(spec.limits))
     stop("specification limits must be provided")
  spec.limits <- as.vector(spec.limits)[1:2]
  LSL <- spec.limits[1]
  if(!(is.numeric(LSL) & is.finite(LSL))) LSL <- NA
  USL <- spec.limits[2]
  if(!(is.numeric(USL) & is.finite(USL))) USL <- NA
  if(is.na(LSL) & is.na(USL))
     stop("invalid specification limits")

  has.target <- (!missing(target))
  if(!has.target) 
    { target <- mean(spec.limits, na.rm=TRUE)
      if(!is.na(LSL) & !is.na(USL)) has.target <- TRUE
  }
     
  if (is.na(LSL))
     { if (target > USL)
           warning("target value larger than one-sided specification limit...") }
  if (is.na(USL))
     { if (target < LSL)
           warning("target value smaller than one-sided specification limit...") }
  if (!is.na(LSL) & !is.na(USL))
     { if (target < LSL || target > USL)
       warning("target value is not within specification limits...") }
       
  if (missing(nsigmas))
     if (is.null(object$nsigmas))
        stop("nsigmas not available in the 'qcc' object. Please provide nsigmas.") 
     else  nsigmas <- object$nsigmas
  
  if (confidence.level < 0 | confidence.level > 1)
     stop("the argument confidence.level must be a value between 0 and 1") 

  # computes process capability indices
  Cp <- (USL - LSL) / (2*nsigmas*std.dev)
  Cp.u <- (USL-center)/(nsigmas*std.dev)
  Cp.l <- (center-LSL)/(nsigmas*std.dev)
  Cp.k <- min(Cp.u, Cp.l)
  # Cpm <- (USL - LSL) / (2*nsigmas*sqrt(sum((x-target)^2)/(n-1)))
  Cpm <- Cp / sqrt(1+((center-target)/std.dev)^2)

  # compute confidence limits 
  alpha <- 1-confidence.level
  Cp.limits   <- Cp*sqrt(qchisq(c(alpha/2, 1-alpha/2), n-1)/(n-1))
  Cp.u.limits <- Cp.u*(1+c(-1,1)*qnorm(confidence.level)*
                       sqrt(1/(9*n*Cp.u^2)+1/(2*(n-1))))
  Cp.l.limits <- Cp.l*(1+c(-1,1)*qnorm(confidence.level)*
                       sqrt(1/(9*n*Cp.l^2)+1/(2*(n-1))))
  Cp.k.limits <- Cp.k*(1+c(-1,1)*qnorm(1-alpha/2)*
                       sqrt(1/(9*n*Cp.k^2)+1/(2*(n-1))))
  df <- n*(1+((center-target)/std.dev)^2)/
          (1+2*((center-target)/std.dev)^2)
  Cpm.limits <- Cpm*sqrt(qchisq(c(alpha/2, 1-alpha/2), df)/df)
  names(Cp.limits) <- names(Cp.k.limits) <- names(Cpm.limits) <- 
    c(paste(round(100*alpha/2, 1), "%", sep=""),
      paste(round(100*(1-alpha/2), 1), "%", sep=""))

  if(is.na(LSL))  exp.LSL <- NA
  else { exp.LSL <- pnorm((LSL-center)/std.dev) * 100
         if(exp.LSL < 0.01) exp.LSL <- 0 }
  if(is.na(USL))  exp.USL <- NA
  else { exp.USL <- (1-pnorm((USL-center)/std.dev)) * 100
         if(exp.USL < 0.01) exp.USL <- 0 }
  obs.LSL <- sum(x<LSL)/n * 100
  obs.USL <- sum(x>USL)/n * 100
  
  tab <- cbind(c(Cp, Cp.l, Cp.u, Cp.k, Cpm),
               rbind(Cp.limits, Cp.l.limits, Cp.u.limits, 
                     Cp.k.limits, Cpm.limits))
  rownames(tab) <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm")
  colnames(tab) <- c("Value", names(Cp.limits))

  out <- list(data = x, data.name = object$data.name,
              center = center, std.dev = std.dev, 
              has.target = has.target, target = target, 
              spec.limits = { sl <- c(LSL, USL)
                              names(sl) <- c("LSL", "USL")
                              sl },
              indices = tab, 
              exp = { exp <- c(exp.LSL, exp.USL)/100
                      names(exp) <- c("Exp < LSL", "Exp > USL")
                      exp }, 
              obs = { obs <- c(obs.LSL, obs.USL)/100
                      names(obs) <- c("Obs < LSL", "Obs > USL")
                      obs } )
  class(out) <- "processCapability"
  if(plot) plot(out, ...) 

  return(out)
}

print.processCapability <- function(x, digits = getOption("digits"), ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Process Capability Analysis"), 
                width = min(getOption("width"),50)), "\n\n")
  
  cat(paste(formatC("Number of obs = ", width=16),
            formatC(length(object$data), width=12, flag="-"),
            formatC("Target = ", width=10),
            ifelse(object$has.target, 
                   formatC(object$target, digits=digits, flag="-"), ""),
            "\n", sep=""))
  cat(paste(formatC("Center        = ", width=16),
            formatC(object$center, digits=digits, width=12, flag="-"),
            formatC("LSL    = ", width=10),
            ifelse(is.na(object$spec.limits[1]), "",
                   formatC(object$spec.limits[1], digits=digits, flag="-")),
            "\n", sep=""))
  cat(paste(formatC("StdDev        = ", width=16),
            formatC(object$std.dev, digits=digits, width=12, flag="-"),
            formatC("USL    = ", width=10),
            ifelse(is.na(object$spec.limits[2]), "",
                   formatC(object$spec.limits[2], digits=digits, flag="-")),
            "\n", sep=""))
            
  indices <- object$indices
  names(dimnames(indices)) <- c("Capability indices", "")
  print(indices, digits = 3, na.print="", print.gap=2)
  
  cat("\n")
  cat(paste("Exp<LSL", ifelse(is.na(object$exp[1]), "\t", 
                              paste(signif(object$exp[1], digits=2), "%\t", sep="")), 
            "Obs<LSL", ifelse(is.na(object$obs[1]), "", 
                              paste(signif(object$obs[1], digits=2), "%\n", sep=""))))
  cat(paste("Exp>USL", ifelse(is.na(object$exp[2]), "\t", 
                              paste(signif(object$exp[2], digits=2), "%\t", sep="")),
            "Obs>USL", ifelse(is.na(object$obs[2]), "", 
                              paste(signif(object$obs[2], digits=2), "%\n", sep=""))))
  
  invisible()
}

summary.processCapability <- function(object, ...) 
  print.processCapability(object, ...)

plot.processCapability <- function(x, 
                                   add.stats = qcc.options("add.stats"),
                                   breaks = "scott", 
                                   col = adjustcolor(qcc.options("zones")$fill, alpha.f = 0.5), 
                                   border = "white",
                                   digits = getOption("digits"),
                                   restore.par = TRUE, ...) 
{
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "processCapability")))
     stop("an object of class `processCapability' is required")

  xlim <- range(object$data, object$spec.limits, object$target, na.rm = TRUE)
  xlim <- extendrange(r = xlim, f = 0.1)
  x  <- seq(min(xlim), max(xlim), length=250)
  dx <- dnorm(x, object$center, object$std.dev)
  h <- hist(object$data, breaks = breaks, plot=FALSE) # compute histogram
  ylim <- range(h$density, dx)
  ylim <- ylim+diff(ylim)*c(0,0.05)
  nobs <- length(object$data)
  Cp   <- object$indices[1,1]
  Cp_l <- object$indices[2,1]
  Cp_u <- object$indices[3,1]
  Cp_k <- object$indices[4,1]
  Cpm  <- object$indices[5,1]
  
  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  cex.labels <- par("cex")*qcc.options("cex")
  cex.stats <- par("cex")*qcc.options("cex.stats")
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mgp = c(2.1, 0.8, 0),
      mar = c(4.1,2.1,1.1,2.1),
      oma = if(add.stats) c(5.5*cex.stats, 0, 1.5*cex.labels, 0) 
            else          c(0, 0, 1.5*cex.labels, 0))

  plot(0, 0, type="n", xlim = xlim, ylim = ylim,
       axes = FALSE, ylab="", xlab = object$data.name)
  usr <- par()$usr
  rect(usr[1], usr[3], usr[2], usr[4], col = qcc.options("bg.figure"))
  axis(1, cex.axis = par("cex.axis")*0.9)
  box()
  mtext("Process capability analysis", 
        side = 3, outer = TRUE, 
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"),
        col  = par("col.main"))
  # draw histogram
  plot(h, add = TRUE, freq = FALSE, 
       col = col, border = border) 
  # add graphical info
  abline(v = object$spec.limits, lty = 3, lwd = 1)
  text(object$spec.limits[1], usr[4], "LSL", 
       col = gray(0.3), pos = 3, offset = 0.2, 
       cex = par("cex")*qcc.options("cex.stats"), xpd = TRUE)
  text(object$spec.limits[2], usr[4], "USL", 
       col = gray(0.3), pos = 3, offset = 0.2, 
       cex = par("cex")*qcc.options("cex.stats"), xpd = TRUE)
  if(object$has.target)
  { 
    abline(v = object$target, lty = 2, lwd = 1)
    text(object$target, usr[4], "Target", 
         col = gray(0.3), pos = 3, offset = 0.2, 
         cex = par("cex")*qcc.options("cex.stats"), xpd = TRUE)
  }
  lines(x, dx, lty=1)

  if(add.stats) 
    { 
      at <- c(0.07, 0.35, 0.56, 0.75)
      # write info at bottom
      #--
      mtext(paste("Number of obs = ", nobs, sep = ""), 
            side = 1, outer = TRUE, line = 0, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Center = ", signif(object$center, digits), sep = ""), 
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("StdDev = ", signif(object$std.dev, digits), sep = ""), 
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      #--
      mtext(ifelse(object$has.target, 
                   paste("Target = ", signif(object$target, digits), sep = ""),
                   paste("Target = ")),
            side = 1, outer = TRUE, line = 0, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("LSL = ", ifelse(is.na(object$spec.limits[1]), "", 
                                   signif(object$spec.limits[1], digits)), sep = ""), 
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("USL = ", ifelse(is.na(object$spec.limits[2]), "", 
                                   signif(object$spec.limits[2], digits)), sep = ""), 
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      #--
      mtext(paste("Cp     = ", ifelse(is.na(Cp), "", signif(Cp, 3)), sep = ""), 
            side = 1, outer = TRUE, line = 0, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cp_l  = ", ifelse(is.na(Cp_l), "", 
                                     signif(Cp_l, 3)), sep = ""), 
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cp_u = ", ifelse(is.na(Cp_u), "", 
                                    signif(Cp_u, 3)), sep = ""), 
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cp_k = ", ifelse(is.na(Cp_k), "", 
                                    signif(Cp_k, 3)), sep = ""), 
            side = 1, outer = TRUE, line = 3*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cpm  = ", ifelse(is.na(Cpm), "", signif(Cpm, 3)), sep = ""), 
            side = 1, outer = TRUE, line = 4*cex.stats, adj = 0, at = at[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      #--
      mtext(paste("Exp<LSL ", ifelse(is.na(object$exp[1]), "", 
                                     paste(signif(object$exp[1], 2), "%", sep="")), sep = ""), 
            side = 1, outer = TRUE, line = 0, adj = 0, at = at[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Exp>USL ", ifelse(is.na(object$exp[2]), "", paste(signif(object$exp[2], 2), "%", sep="")), sep = ""),
            side = 1, outer = TRUE, line = 1*cex.stats, adj = 0, at = at[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Obs<LSL ", ifelse(is.na(object$obs[1]), "", 
                                     paste(signif(object$obs[1], 2), "%", sep="")), sep = ""), 
            side = 1, outer = TRUE, line = 2*cex.stats, adj = 0, at = at[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Obs>USL ", ifelse(is.na(object$obs[2]), "", 
                                     paste(signif(object$obs[2], 2), "%", sep="")), sep = ""), 
            side = 1, outer = TRUE, line = 3*cex.stats, adj = 0, at = at[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
  }

  invisible()
}
