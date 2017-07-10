#-------------------------------------------------------------------#
#                                                                   #
#                     PARETO CHART                                  #
#                                                                   #
#-------------------------------------------------------------------#

pareto.chart <- function(data, plot = TRUE, ...)
{ 
  call <- match.call(expand.dots = TRUE)
  varname <- deparse(substitute(data))
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
  names(dimnames(tab)) <- c("", paste("\nPareto chart analysis for", varname))
  tab <- as.table(tab)
  class(tab) <- c("pareto.chart", class(tab))
  attr(tab, "call") <- call
  attr(tab, "varname") <- varname
  if(plot) plot(tab, ...)
  return(tab)
}

print.pareto.chart <- function(x, ...) print.table(x, ...)

plot.pareto.chart <- function(x, xlab = NULL, 
                              ylab = "Frequency", 
                              ylab2 = "Cumulative Percentage", 
                              cumperc = seq(0, 100, by = 25), 
                              ylim = NULL, 
                              main = NULL, 
                              col = blues.colors(nlevels), 
                              ...)
{
  call <- attr(x, "call")
  nlevels <- nrow(x)
  freq <- x[,1]
  cumfreq <- x[,2]
  cumperc <- cumperc[cumperc >= 0 & cumperc <= 100]
  q <- quantile(seq(0, max(cumfreq, na.rm = TRUE), 
                    by = max(cumfreq, na.rm = TRUE) / 100), 
                cumperc/100)
  if(is.null(ylim)) ylim <- c(0, max(cumfreq, na.rm = TRUE)*1.05)
  if(is.null(main)) main <- paste("Pareto Chart for", attr(x, "varname"))
  # set las and mar if not provided by user
  w <- max(sapply(rownames(x), nchar))
  if(is.null(call$las)) las <- 3 else las <- call$las
  if(is.null(call$mar))
    { if (las==1) mar <- c(1,1,0,2)  
      else        mar <- c(log(max(w),2),0,0,2) }
  else mar <- call$mar
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"),
      mar = pmax(par("mar")+mar,c(4.1,4.1,3.1,4.1)), 
      las = las, 
      cex = oldpar$cex*qcc.options("cex"))
  #
  pc <- barplot(x[,1], width = 1, space = 0.2, col = col,
                ylim = ylim, ylab = ylab, xlab = xlab, yaxt = "n",  
                ...)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  box()
  top.line <- par("mar")[3]/3
  mtext(main, side = 3, line = top.line, las = 1,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  # adding line for percentage level overwrite bars...
  abline(h = q, col = "lightgrey", lty = 3)
  # ... so we redraw bars (not nice but works!)
  rect(pc-0.5, rep(0,nlevels), pc+0.5, freq, col = col)
  lines(pc, cumfreq, type = "b", cex = 0.8*par("cex"), pch = 19)
  box()
  axis(2, las = 3)
  axis(4, at = q, las = 3, labels = paste(cumperc, "%", sep = ""))
  mtext(ylab2, side = 4, line = 2.5, las = 3, 
        cex = par("cex")*qcc.options("cex.stats"))
}


