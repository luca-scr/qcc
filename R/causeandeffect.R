#-------------------------------------------------------------------#
#                                                                   #
#                  Cause-and-Effect diagram                         #
#                                                                   #
#-------------------------------------------------------------------#

cause.and.effect <- function(cause, effect, title = "Cause-and-Effect diagram", cex=c(1,0.9,1), font=c(1,3,2))
{

  # running mean of successive pairs of obs
  mean2 <- function(x)
  { m <- rep(NA, length(x)-1)
    for (i in 1:(length(x)-1))
         m[i] <- mean(x[c(i,i+1)])
    return(m)
  }

  nc <- length(cause)
  ncup <- nc - round(nc/2)
  nclo <- nc - ncup
  ncc <- max(sapply(cause, length))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(bg  = qcc.options("bg.margin"), 
      cex = oldpar$cex * qcc.options("cex"),
      mar = c(2,2,3,2))

  plot(0:100, 0:100, type = "n", xlab = "", ylab = "", axes = FALSE)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
       col = qcc.options("bg.figure"))
  # box()
  mtext(title, side = 3, line = par("mar")[3]/3,
        font = par("font.main"), 
        cex  = qcc.options("cex"), 
        col  = par("col.main"))
  
  usr <- par("usr")
  we <- strwidth(effect, units="user")*1.1
  wc <- max(unlist(sapply(cause, strwidth, units="user")))
  hc <- max(strheight(effect, units="user"),
            unlist(sapply(cause, strheight, units="user")))

  # draw effect
  arrows(0, 50, usr[2]-we-1, 50, code=2, length=0.1, angle=20)
  text(usr[2]-we,50, effect, adj=c(0,0.5), cex=cex[3], font=font[3])

  # draw branches and cause labels
  a <- (usr[2]-we)/(max(ncup,nclo)+1)
  ac <-  a*(0:(max(ncup,nclo)))
  for (i in 1:(length(ac)-1))
      { segments(mean2(ac)[i], 95, ac[i+1], 50) 
        text(mean2(ac)[i], 96, names(cause)[i], 
             pos=3, offset=0.5, cex=cex[1], font=font[1]) 
        if (i <= nclo)
           { segments(mean2(ac)[i], 5, ac[i+1], 50)
             text(mean2(ac)[i], 4, names(cause)[[ncup+i]], 
                  pos=1, offset=0.5, cex=cex[1], font=font[1]) }
      }

  # draw labels for upper branches
  for (j in 1:ncup)
      { b <- (50-95)/(ac[j+1]-mean2(ac)[j])
        a <- 95-b*mean2(ac)[j]
        y <- rev(50+cumsum((95-50)/(ncc+1))*(1:(ncc)))
        x <- (y-a)/b
        for (i in 1:length(y))
            { label <- cause[[j]][i]
              if (!is.na(label))
                 text(x[i], y[i], label, pos=4, 
                      offset=0.2, cex=cex[2], font=font[2]) 
             }
      }       
  # draw labels for lower branches
  for (j in 1:ncup)
      { b <- (50-5)/(ac[j+1]-mean2(ac)[j])
        a <- 5-b*mean2(ac)[j]
        y <- cumsum((95-50)/(ncc+1))*(1:(ncc))
        x <- (y-a)/b
        if (j <= nclo)
           for (i in 1:length(y))
               { label <- cause[[ncup+j]][i]
                 if (!is.na(label))
                    text(x[i], y[i], label, pos=4, 
                         offset=0.2, cex=cex[2], font=font[2])
               }
      }

  invisible()
}
