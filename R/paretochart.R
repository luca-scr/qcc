#-------------------------------------------------------------------#
#                                                                   #
#                     PARETO CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

paretoChart <- function(data, plot = TRUE, ...)
{ 
  call <- match.call(expand.dots = TRUE)
  data.name <- deparse(substitute(data))
  data <- as.table(data)
  if(length(dim(data))>1) 
    stop("only one-dimensional object (table, vector, etc.) may be provided")
  # 
  data <- sort(data, decreasing = TRUE, na.last = TRUE)
  csum <- cumsum(data)
  tab <- cbind(data, csum, 
               data/max(csum, na.rm = TRUE)*100, 
               csum/max(csum, na.rm = TRUE)*100) 
  colnames(tab) <- c("Frequency", "Cum.Freq.", 
                     "Percentage", "Cum.Percent.")
  names(dimnames(tab)) <- c(data.name, "")
  
  # create object of class 'paretoChart'
  object <- list(call = call, 
                 data.name = data.name, 
                 tab = tab)
  class(object) <- "paretoChart"
  
  if(plot) plot(object, ...)
  return(object)
}

print.paretoChart <- function(x, digits = getOption("digits") - 3, ...)
{
  object <- x   # Argh.  Really want to use 'object' anyway
  cat(cli::rule(left = crayon::bold("Pareto Chart"), 
                width = min(getOption("width"),50)), "\n")
  print(object$tab, digits = digits, ...)
}

plot.paretoChart <- function(x, xlab = NULL, 
                             ylab = "Frequency", 
                             ylab2 = "Cumulative Percentage", 
                             cumperc = seq(0, 100, by = 25), 
                             ylim = NULL, 
                             main = NULL, 
                             col = blues.colors(nlevels), 
                             ...)
{
  call <- x$call
  freq <- x$tab[,1]
  cumfreq <- x$tab[,2]
  cumperc <- cumperc[cumperc >= 0 & cumperc <= 100]
  nlevels <- length(freq)
  q <- quantile(seq(0, max(cumfreq, na.rm = TRUE), 
                    by = max(cumfreq, na.rm = TRUE) / 100), 
                cumperc/100)
  if(is.null(ylim)) ylim <- c(0, max(cumfreq, na.rm = TRUE)*1.05)
  if(is.null(main)) main <- paste("Pareto Chart for", x$data.name)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  # plot.new()

  # set las and mar if not provided by user
  # browser()
  w <- max(strwidth(rownames(x$tab), units = "inch"))
  las <- if(is.null(call$las)) 3 else eval(call$las)
  if(is.null(call$mar))
  { 
    c <- max(par("mar")/par("mai"))
    mar <- c(3.1,4.1,1,4.1)
    if(las==2 | las==3) mar[1] <- c*w+1
  } else 
  { 
    mar <- eval(call$mar)
  }
  if(!is.null(xlab)) mar[1] <- mar[1] + 1.5
  cex.labels <- par("cex")*qcc.options("cex")

  par(bg  = qcc.options("bg.margin"),
      oma = c(0, 0, 1.5*cex.labels, 0),
      mar = mar) 
  #
  pc <- barplot(freq, width = 1, space = 0.2, col = col,
                ylim = ylim, ylab = ylab, xlab = xlab, yaxt = "n", 
                cex.names = par("cex.axis")*0.9, 
                cex.axis = par("cex.axis")*0.9, 
                cex.lab = par("cex.axis")*0.9, 
                las = las, ...)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  
  # box()
  mtext(main, side = 3, outer = TRUE, 
        line = 0, adj = 0, at = par("plt")[1],
        font = par("font.main"), 
        cex  = par("cex")*qcc.options("cex"),
        col  = par("col.main"))
  # adding line for percentage level overwrite bars...
  abline(h = q, col = "lightgrey", lty = 3)
  # ... so we redraw bars (not nice but works!)
  rect(pc-0.5, rep(0,nlevels), pc+0.5, freq, col = col)
  lines(pc, cumfreq, type = "b", cex = 0.8*par("cex"), pch = 19)
  box()
  axis(2, las = las, cex.axis = par("cex.axis")*0.9)
  axis(4, at = q, las = las, labels = paste(cumperc, "%", sep = ""),
       cex.axis = par("cex.axis")*0.9)
  mtext(ylab2, side = 4, line = 2.5, las = las, 
        cex = par("cex.axis")*0.9)
}
