#-------------------------------------------------------------------#
#                                                                   #
#    Process Capability Analysis                                    #
#                                                                   #
#-------------------------------------------------------------------#

process.capability <- function(object, spec.limits, target, std.dev, nsigmas, confidence.level = 0.95, breaks="scott", add.stats=TRUE, print=TRUE, digits = getOption("digits"), restore.par=TRUE)
{
# Computes process capability indices for a qcc object of type "xbar" 
# and plot the histogram

  if ((missing(object)) | (!inherits(object, "qcc")))
     stop("an object of class 'qcc' is required")
  if (!(object$type=="xbar" | object$type=="xbar.one"))
     stop("Process Capability Analysis only available for charts type \"xbar\" and \"xbar.one\" charts")

  x <- as.vector(object$data)
  x <- x[!is.na(x)]
  sizes <- object$sizes
  center <- object$center
  if (missing(std.dev))
     std.dev <- object$std.dev
  n <- length(x)
  title <- paste("Process Capability Analysis\nfor", object$data.name)

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
  xlim <- range(x, USL, LSL, target, na.rm = TRUE)
  xlim <- xlim+diff(xlim)*c(-0.1,0.1)
  xx <- seq(min(xlim), max(xlim), length=250)
  dx <- dnorm(xx, center, std.dev)
  h <- hist(x, breaks = breaks, plot=FALSE) # compute histogram
  ylim <- range(h$density, dx)
  ylim <- ylim+diff(ylim)*c(0,0.05)

  tab <- cbind(c(Cp, Cp.l, Cp.u, Cp.k, Cpm),
               rbind(Cp.limits, Cp.l.limits, Cp.u.limits, 
                     Cp.k.limits, Cpm.limits))
  rownames(tab) <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm")
  colnames(tab) <- c("Value", names(Cp.limits))

  oldpar <- par(no.readonly = TRUE)
  if(restore.par) on.exit(par(oldpar))
  mar <- c(4.1,2.1,3.6,2.1)
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = if(add.stats) pmax(mar, c(8.6+is.null(center)*-1,0,0,0)) else mar)

  plot(0, 0, type="n", xlim = xlim, ylim = ylim,
       axes = FALSE, ylab="", xlab = "")
  usr <- par()$usr
  rect(usr[1], usr[3], usr[2], usr[4], col = qcc.options("bg.figure"))
  axis(1); box()
  top.line <- par("mar")[3]-length(capture.output(cat(title))) - 0.5
  mtext(title, side = 3, line = top.line,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  # draw histogram
  plot(h, add = TRUE, freq = FALSE) 
  # add graphical info
  abline(v=c(LSL,USL), col=2, lty=3, lwd=2)
  text(LSL, usr[4], "LSL", pos=3, offset=0.2, cex=0.8, xpd = TRUE)
  text(USL, usr[4], "USL", pos=3, offset=0.2, cex=0.8, xpd = TRUE)
  if(has.target)
    { abline(v=target, col=2, lty=2, lwd=2)
      text(target, usr[4], "Target", pos=3, offset=0.2, cex=0.8, xpd = TRUE) }
  lines(xx, dx, lty=2)

  if(add.stats) 
    { # computes the x margins of the figure region
      plt <- par()$plt
      px <- diff(usr[1:2])/diff(plt[1:2])
      xfig <- c(usr[1]-px*plt[1], usr[2]+px*(1-plt[2]))
      at.col <- xfig[1] + diff(xfig[1:2])*c(0.07, 0.35, 0.56, 0.75)
      top.line <- 3
      # write info at bottom
      #--
      mtext(paste("Number of obs = ", n, sep = ""), 
            side = 1, line = top.line, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Center = ", signif(center, digits), sep = ""), 
            side = 1, line = top.line+1, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("StdDev = ", signif(std.dev, digits), sep = ""), 
            side = 1, line = top.line+2, adj = 0, at = at.col[1],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      #--
      mtext(ifelse(has.target, paste("Target = ", signif(target, digits), sep = ""),
                               paste("Target = ")),
            side = 1, line = top.line, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("LSL = ", ifelse(is.na(LSL), "", signif(LSL, digits)), sep = ""), 
            side = 1, line = top.line+1, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("USL = ", ifelse(is.na(USL), "", signif(USL, digits)), sep = ""), 
            side = 1, line = top.line+2, adj = 0, at = at.col[2],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      #--
      mtext(paste("Cp     = ", ifelse(is.na(Cp), "", signif(Cp, 3)), sep = ""), 
            side = 1, line = top.line, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cp_l  = ", ifelse(is.na(Cp.l), "", signif(Cp.l, 3)), sep = ""), 
            side = 1, line = top.line+1, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cp_u = ", ifelse(is.na(Cp.u), "", signif(Cp.u, 3)), sep = ""), 
            side = 1, line = top.line+2, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cp_k = ", ifelse(is.na(Cp.k), "", signif(Cp.k, 3)), sep = ""), 
            side = 1, line = top.line+3, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Cpm  = ", ifelse(is.na(Cpm), "", signif(Cpm, 3)), sep = ""), 
            side = 1, line = top.line+4, adj = 0, at = at.col[3],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      #--
      mtext(paste("Exp<LSL ", ifelse(is.na(exp.LSL), "", paste(signif(exp.LSL, 2), "%", sep="")), sep = ""), 
            side = 1, line = top.line, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Exp>USL ", ifelse(is.na(exp.USL), "", paste(signif(exp.USL, 2), "%", sep="")), sep = ""),
            side = 1, line = top.line+1, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Obs<LSL ", ifelse(is.na(obs.LSL), "", paste(signif(obs.LSL, 2), "%", sep="")), sep = ""), 
            side = 1, line = top.line+2, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
      mtext(paste("Obs>USL ", ifelse(is.na(obs.USL), "", paste(signif(obs.USL, 2), "%", sep="")), sep = ""), 
            side = 1, line = top.line+3, adj = 0, at = at.col[4],
            font = qcc.options("font.stats"),
            cex = par("cex")*qcc.options("cex.stats"))
    }

  if(print)
    { cat("\nProcess Capability Analysis\n")
      cat("\nCall:\n", deparse(match.call()), "\n\n", sep = "")
      cat(paste(formatC("Number of obs = ", width=16), 
                formatC(n, width=12, flag="-"), 
                formatC("Target = ", width=10), 
                ifelse(has.target, formatC(signif(target,digits=digits), flag="-"), ""),
                "\n", sep=""))
      cat(paste(formatC("Center = ", width=16), 
                formatC(signif(center, digits=digits), width=12, flag="-"),
                formatC("LSL = ", width=10), 
                ifelse(is.na(LSL), "", formatC(signif(LSL, digits=digits), flag="-")),
                "\n", sep=""))
      cat(paste(formatC("StdDev = ", width=16), 
                formatC(signif(std.dev, digits=digits), width=12, flag="-"),
                formatC("USL = ", width=10), 
                ifelse(is.na(USL), "", formatC(signif(USL, digits=digits), flag="-")),
                "\n", sep=""))
      cat("\nCapability indices:\n\n")
      print(tab, digits=4, na.print="", print.gap=2)
      cat("\n")
      cat(paste("Exp<LSL", ifelse(is.na(exp.LSL), "\t", 
                                   paste(format(exp.LSL, digits=2), "%\t", sep="")), 
                "Obs<LSL", ifelse(is.na(obs.LSL), "", 
                                   paste(format(obs.LSL, digits=2), "%\n", sep=""))))
      cat(paste("Exp>USL", ifelse(is.na(exp.USL), "\t", 
                                   paste(format(exp.USL, digits=2), "%\t", sep="")),
                "Obs>USL", ifelse(is.na(obs.USL), "", 
                                   paste(format(obs.USL, digits=2), "%\n", sep=""))))
    }

  invisible(list(nobs = n, center = center, std.dev = std.dev, 
                 target = target, 
                 spec.limits = { sl <- c(LSL, USL)
                                 names(sl) <- c("LSL", "USL")
                                 sl },
                 indices = tab, 
                 exp = { exp <- c(exp.LSL, exp.USL)/100
                         names(exp) <- c("Exp < LSL", "Exp > USL")
                         exp }, 
                 obs = { obs <- c(obs.LSL, obs.USL)/100
                         names(obs) <- c("Obs < LSL", "Obs > USL")
                         obs }
                 ))
}

