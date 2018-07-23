.ddPlot.MatNtNtCoCo <- function(data, ...,
                                dist.x = NormType(), 
                                dist.y  = NormType(),
                                cutoff.x = cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x),
                                cutoff.y = cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y),
                                cutoff.quantile.x = 0.95,  
                                cutoff.quantile.y = cutoff.quantile.x,
                                transform.x,
                                transform.y = transform.x,
                                id.n,
                                cex.pts = 1,
                                lab.pts,
                                jitter.pts = 0,
                                alpha.trsp = NA,
                                adj =0,
                                cex.idn = 1,
                                col.idn = par("col"),
                                lty.cutoff,
                                lwd.cutoff,
                                col.cutoff = "red",
                                text.abline = TRUE,
                                text.abline.x = NULL, text.abline.y = NULL,
                                cex.abline = par("cex"), col.abline = col.cutoff,
                                font.abline = par("font"), adj.abline = c(0,0),
                                text.abline.x.x = NULL, text.abline.x.y = NULL,
                                text.abline.y.x = NULL, text.abline.y.y = NULL,
                                text.abline.x.fmt.cx = "%7.2f",
                                text.abline.x.fmt.qx = "%4.2f%%",
                                text.abline.y.fmt.cy = "%7.2f",
                                text.abline.y.fmt.qy = "%4.2f%%",
                                jitter.fac = 10,
                                jitter.tol = .Machine$double.eps,
                                doplot = TRUE){

       mc <- match.call(expand.dots = FALSE)
       dots <- mc$"..."

       if(missing(jitter.pts)||is.null(jitter.pts)) jitter.pts <- 0
       jitter.pts <- rep(jitter.pts,length.out=2)
       if(missing(jitter.tol)||is.null(jitter.tol)) jitter.tol <- .Machine$double.eps
       jitter.tol <- rep(jitter.tol,length.out=2)
       if(missing(jitter.fac)||is.null(jitter.fac)) jitter.fac <- 10
       jitter.fac <- rep(jitter.fac,length.out=2)

       col <- if(is.null(dots$col)) par("col") else dots$col
       if(!is.na(alpha.trsp)) col <- addAlphTrsp2col(col, alpha.trsp)


       id.n1 <- 1:ncol(data)

       if(missing(id.n) || is.null(id.n))
          id.n <- id.n1


       if(missing(lab.pts)|| is.null(lab.pts)){
          lab.pts <-  if(!is.null(colnames(data))) colnames(data) else id.n1
       }

       data <- data[,id.n, drop = FALSE]
       lab.pts <- lab.pts[id.n]
       id.n1 <- id.n1[id.n]

       id.n0 <- 1:length(id.n)

       data.x <- if(missing(transform.x)||is.null(transform.x)) data else
                    transform.x(data)
                    #if(!is(IC.x,"IC")) stop("Argument 'IC.x' of 'ddPlot' must be an 'IC'")
                    #   apply(data,2,function(xx) evalIC(IC.x,xx))}
       data.y <- if(missing(transform.y)||is.null(transform.y)) data else
                    transform.y(data)
                    #if(!is(IC.y,"IC")) stop("Argument 'IC.y' of 'ddPlot' must be an 'IC'")
                    #   apply(data,2,function(xx) evalIC(IC.y,xx))}


      if(is.null(dist.x)) dist.x <- NormType()
      if(is.null(dist.y)) dist.y <- NormType()

      if(is.null(dots$xlab)) dots$xlab <- name(dist.x)
      if(is.null(dots$ylab)) dots$ylab <- name(dist.y)

      if(!is.null(dots$log)){
	          if(grepl("x",dots$log)) dots$xlab <- paste(dots$xlab, "(log-scale)",
	                                               sep="  ")
	          if(grepl("y",dots$log)) dots$ylab <- paste(dots$ylab, "(log-scale)",
	                                               sep="  ")
	       }

      if(is.null(cutoff.quantile.x))
         cutoff.quantile.x <- 0.95

      if(is.null(cutoff.quantile.y))
         cutoff.quantile.y <- cutoff.quantile.x

      if(is.null(cutoff.x))
         cutoff.x <- cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x)
      else {assign("norm", dist.x, envir=environment(fct(cutoff.x)))
            assign("cutoff.quantile", cutoff.quantile.x, envir=environment(fct(cutoff.x)))
            assign("..trf", if(missing(transform.x)||is.null(transform.x)) function(x)x else transform.x,
                   envir=environment(fct(cutoff.x)))}

      if(is.null(cutoff.y))
         cutoff.y <- cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y)
      else {assign("norm", dist.y, envir=environment(fct(cutoff.y)))
            assign("cutoff.quantile", cutoff.quantile.y, envir=environment(fct(cutoff.y)))
            assign("..trf", if(missing(transform.y)||is.null(transform.y)) function(x)x else transform.y,
                   envir=environment(fct(cutoff.y)))}

      if(!is(dist.x, "NormType")) stop("Argument 'dist.x' of 'ddPlot' must be of class 'NormType'")
      if(!is(dist.y, "NormType")) stop("Argument 'dist.y' of 'ddPlot' must be of class 'NormType'")
      if(!is(cutoff.x, "cutoff")) stop("Argument 'cutoff.x' of 'ddPlot' must be of class 'cutoff'")
      if(!is(cutoff.y, "cutoff")) stop("Argument 'cutoff.y' of 'ddPlot' must be of class 'cutoff'")

      ndata.x <- fct(dist.x)(data.x)
      ndata.y <- fct(dist.y)(data.y)
      
#      print(head(ndata.x))

      co.x <- fct(cutoff.x)(data.x)
      co.y <- fct(cutoff.y)(data.y)

      if(is.null(adj)) adj <- 0
      if(missing(cex.idn)||is.null(cex.idn))
         cex.idn <- if(is.null(dots$cex)) 1 else dots$cex

      if(missing(col.idn)||is.null(col.idn))
         col.idn <- if(is.null(dots$col)) par("col") else dots$col

      if(is.null(dots$lwd)) dots$lwd <- par("lwd")
      if(is.null(dots$lty)) dots$lty <- par("lty")

      if(missing(col.cutoff) || is.null(col.cutoff)) col.cutoff <- "red"
      col.cutoff <- rep(col.cutoff,length.out=2)
      if((missing(lty.cutoff)|| is.null(lty.cutoff)) && !is.null(dots$lty)) lty.cutoff <- dots$lty
      if((missing(lwd.cutoff)|| is.null(lwd.cutoff)) && !is.null(dots$lwd)) lwd.cutoff <- dots$lwd
      if((missing(cex.abline)|| is.null(cex.abline)) && !is.null(dots$cex)) cex.abline <- dots$cex
      if((missing(cex.abline)|| is.null(cex.abline))) cex.abline <- par("cex")
      if((missing(adj.abline)|| is.null(adj.abline)) && !is.null(dots$adj)) adj.abline <- dots$adj
      if((missing(adj.abline)|| is.null(adj.abline))) adj.abline <- c(0.5,0.5)
      if((missing(font.abline)|| is.null(font.abline)) && !is.null(dots$font)) font.abline <- dots$font
      if((missing(font.abline)|| is.null(font.abline))) font.abline <- par("font")

      pdots <- .makedotsLowLevel(dots)
      pdots$pch <- if(is.null(dots$pch)) "." else dots$pch
      pdots$cex <- cex.pts
      pdots$xlab <- dots$xlab
      pdots$ylab <- dots$ylab
      pdots$nsim <- NULL
      pdots$x <- NULL
      pdots$y <- NULL
      pdots$offset <- NULL
      pdots$pos <- NULL
      pdots$untf <- NULL

      abdots <- .makedotsAB(dots)
      if(!missing(lty.cutoff)) abdots$lty <- lty.cutoff[[1]]
      if(!missing(lwd.cutoff)) abdots$lwd <- lwd.cutoff[1]
      if(!missing(col.cutoff)) abdots$col <- col.cutoff[1]

      abdots <- list(abdots,abdots)

      if(!is.null(abdots$lty))
	          if(is.list(lty.cutoff)) abdots[[2]]$lty <-  lty.cutoff[[2]]
      if(!is.null(abdots$lwd))
	         if(length(lwd.cutoff)>1) abdots[[2]]$lwd <-  lwd.cutoff[2]
      if(!is.null(abdots$col))
	         if(length(col.cutoff)>1) abdots[[2]]$col <-  col.cutoff[2]

      if(missing(text.abline)||is.null(text.abline)) text.abline <- TRUE
      ab.textL <- rep(text.abline,length.out=2)

      abtdots.x <- abtdots.y <- vector("list",0)
	    cex.abline <- rep(cex.abline, length.out = 2)
	    col.abline <- rep(if(!is.null(col.abline))
                          col.abline else "red", length.out = 2)
      font.abline <- rep(font.abline, length.out = 2)
      adj.abline <- matrix(rep(adj.abline,length.out=4),2,2)


      if(is.null(text.abline.x.fmt.cx))  text.abline.x.fmt.cx <- "%7.2f"
      if(is.null(text.abline.x.fmt.qx))  text.abline.x.fmt.qx <- "%4.2f%%"
      if(is.null(text.abline.y.fmt.cy))  text.abline.y.fmt.cy <- "%7.2f"
      if(is.null(text.abline.y.fmt.qy))  text.abline.y.fmt.qy <- "%4.2f%%"
	    .mpresubs <- function(inx)
                    .presubs(inx, c("%qx", "%qy", "%cx", "%cy"),
                        c(gettextf(text.abline.x.fmt.qx,
                             round(cutoff.quantile.x*100,1)),
                          gettextf(text.abline.y.fmt.qy,
                             round(cutoff.quantile.y*100,1)),
                          gettextf(text.abline.x.fmt.cx,
                             round(co.x,2)),
                          gettextf(text.abline.y.fmt.cy,
                          round(co.y,2))))
      if(!is.null(text.abline.x)){abtdots.x$labels <- .mpresubs(text.abline.x)
         }else{
         abtdots.x$labels <- .mpresubs(gettextf("%%qx-cutoff =%%cx"))
         }
      abtdots.x$cex <- cex.abline[1]
	    abtdots.x$col <- col.abline[1]
	    abtdots.x$font <- font.abline[1]
	    abtdots.x$srt <- NULL
	    abtdots.x$adj <- adj.abline[,1]

      abtdots.y$labels <- if(! is.null(text.abline.y))
                       .mpresubs(text.abline.y) else .mpresubs(gettextf(
                              "%%qy-cutoff =%%cy"))
	    abtdots.y$cex <- cex.abline[2]
	    abtdots.y$col <- col.abline[2]
	    abtdots.y$font <- font.abline[2]
	    abtdots.y$srt <- NULL
	    abtdots.y$adj <- adj.abline[,2]

      tdots <- .makedotsT(dots)
      tdots$cex <- cex.idn
      tdots$col <- col.idn
      tdots$offset <- dots$offset
      tdots$pos <- dots$pos
      tdots$adj <- adj

      pdots$log <- dots$log
      pdots$adj <- par("adj")

      adots <- pdots
      adots$col <- pdots$col.axis
      adots$lty <- pdots$lty.axis
      adots$adj <- par("adj")

      pdots$axes <- FALSE
      pdots$adj <- par("adj")
      ####

#      print(quantile(ndata.x))
#      print(co.x)
#      print(fct(cutoff.x))
#      print(dist.x)
#      print(get("norm", environment(fct(cutoff.x))))
#      print(quantile(ndata.y))
#      print(co.y)
#      print(fct(cutoff.y))
#      print(get("norm", environment(fct(cutoff.y))))
##

      if(!is.null(dots$xlim))
          if(is.logical(dots$xlim))
             if(dots$xlim) pdots$xlim <- c(min(ndata.x)/1.03,max(ndata.x,co.x)*1.03)
       if(!is.null(dots$ylim))
          if(is.logical(dots$ylim))
             if(dots$ylim) pdots$ylim <- c(min(ndata.y)/1.03,max(ndata.y,co.y)*1.03)


      id.x <- id.n0[ndata.x >= co.x*.999]
      id.y <- id.n0[ndata.y >= co.y*.999]
      id.xy <- intersect(id.x,id.y)


      id0.xy <- id.n1[id.xy]
      id0.x <- id.n1[id.x]
      id0.y <- id.n1[id.y]

      ndata.x0 <- ndata.x
      ndata.y0 <- ndata.y
      isna <- is.na(ndata.x0)|is.na(ndata.y0)

      if(any(duplicated(ndata.x0[!isna]/jitter.tol[1])))
          ndata.x0[!isna] <- jitter(ndata.x0[!isna], factor=jitter.pts[1])
      if(any(duplicated(ndata.y0[!isna]/jitter.tol[2])))
          ndata.y0[!isna] <- jitter(ndata.y0[!isna], factor=jitter.pts[2])

      pdots$col <- col
      retV <- list(id.x=id0.x, id.y= id0.y, id.xy = id0.xy,
             qtx = quantile(ndata.x), qty = quantile(ndata.y),
             cutoff.x.v = co.x, cutoff.y.v = co.y)

      if(doplot){
        plotInfo<- list("plotArgs"=NULL)

        plotInfo$PlotArgs <- c(list(x = ndata.x0, y=ndata.y0, type = "p"), pdots)
        plotInfo$BoxArgs <- c(adots)

        do.call(plot, args = plotInfo$PlotArgs)
        do.call(box,args=plotInfo$BoxArgs)

        pusr <- par("usr")
        plotInfo$usr <- pusr

        mid.x <- mean(pusr[c(1,2)])
        mid.y <- mean(pusr[c(3,4)])
        abtdots.y$x <- if(is.null(text.abline.y.x)) mid.x else text.abline.y.x
        abtdots.x$y <- if(is.null(text.abline.x.y)) mid.y else text.abline.x.y

        plotInfo$ablineV <- c(list(v=co.x), abdots[[1]])
        plotInfo$ablineH <- c(list(h=co.y), abdots[[2]])
        do.call(abline, args = plotInfo$ablineV)
	      do.call(abline, args = plotInfo$ablineH)

        if(ab.textL[1]){
           plotInfo$abtextV <- c(list(y=co.y*1.03), abtdots.y)
           do.call(text, args = plotInfo$abtextV)
#         do.call(text, args = c(list(co.x-5,mid.y,paste(cutoff.quantile.y*100,"%-cutoff = ",round(co.x,digits=2)),srt=90)))
        }
        if(ab.textL[2]){
           plotInfo$abtextH <- c(list(x=co.x*1.03), abtdots.x,srt=90)
           do.call(text, args = plotInfo$abtextH)
#      do.call(text, args = c(list(mid.x,co.y+5,paste(cutoff.quantile.x*100," %-cutoff = ",round(co.y,digits=2)))))
        }
        if(length(id.xy)){
           plotInfo$Lab <- c(list(jitter(ndata.x[id.xy],factor=jitter.fac[1]),
                                     jitter(ndata.y[id.xy],factor=jitter.fac[2]),
                                labels=lab.pts[id.xy]), tdots)
           do.call(text, args = plotInfo$Lab)
        }
        plotInfo$retV <- retV
        class(plotInfo) <- c("plotInfo","DiagnInfo")
        return(invisible(plotInfo))
          #axis(side=4)
#      axis(side=1)

      }
      return(invisible(retV))
}
