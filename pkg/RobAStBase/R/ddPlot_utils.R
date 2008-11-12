.ddPlot.MatNtNtCoCo <- function(data, ...,  dist.x = NormType(), dist.y  = NormType(),
       cutoff.x = cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x),
       cutoff.y = cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y),
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj =0, cex.idn = 1,
       col.idn = par("col"), lty.cutoff,
       lwd.cutoff, col.cutoff = "red"){

       dots <- match.call(expand.dots = FALSE)$"..."

       id.n1 <- 1:ncol(data)

       if(missing(id.n) || is.null(id.n))
          id.n <- id.n1


       if(missing(lab.pts)|| is.null(lab.pts)){
          lab.pts <-  if(!is.null(colnames(data))) colnames(data) else 1:ncol(data)
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

      if(is.null(cutoff.quantile.x))
         cutoff.quantile.x <- 0.95

      if(is.null(cutoff.quantile.y))
         cutoff.quantile.y <- cutoff.quantile.x

      if(is.null(cutoff.x))
         cutoff.x <- cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x)
      else {assign("norm", dist.x, environment(fct(cutoff.x)))
            assign("cutoff.quantile", cutoff.quantile.x, environment(fct(cutoff.x)))}

      if(is.null(cutoff.y))
         cutoff.y <- cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y)
      else {assign("norm", dist.y, environment(fct(cutoff.y)))
            assign("cutoff.quantile", cutoff.quantile.y, environment(fct(cutoff.y)))}

      if(!is(dist.x, "NormType")) stop("Argument 'dist.x' of 'ddPlot' must be of class 'NormType'")
      if(!is(dist.y, "NormType")) stop("Argument 'dist.y' of 'ddPlot' must be of class 'NormType'")
      if(!is(cutoff.x, "cutoff")) stop("Argument 'cutoff.x' of 'ddPlot' must be of class 'cutoff'")
      if(!is(cutoff.y, "cutoff")) stop("Argument 'cutoff.y' of 'ddPlot' must be of class 'cutoff'")

      ndata.x <- fct(dist.x)(data.x)
      ndata.y <- fct(dist.y)(data.y)


      if(is.null(adj)) adj <- 0
      if(is.null(cex.idn)) cex.idn <- 1
      if(is.null(col.idn)) col.idn <- par("col")
      if(is.null(col.cutoff)) col.cutoff <- "red"

      if(is.null(dots$lwd)) dots$lwd <- par("lwd")
      if(is.null(dots$lty)) dots$lty <- par("lty")


      pdots <- dots
      pdots$type <- NULL
      pdots$x <- NULL
      pdots$y <- NULL
      pdots$offset <- NULL
      pdots$pos <- NULL
      pdots$log <- NULL
      pdots$untf <- NULL

      abdots <- pdots
      abdots$col <- col.cutoff
      if(!missing(lwd.cutoff)) abdots$lwd <- lwd.cutoff
      if(!missing(lty.cutoff)) abdots$lty <- lty.cutoff
      abdots$pos <- NULL
      abdots$untf <- dots$untf
      abdots$adj <- NULL

      adots <- pdots
      adots$col <- pdots$col.axis
      adots$lty <- pdots$lty.axis
      adots$adj <- par("adj")

      tdots <- pdots
      tdots$cex <- cex.idn
      tdots$col <- col.idn
      tdots$offset <- dots$offset
      tdots$pos <- dots$pos
      tdots$adj <- adj

      pdots$axes <- FALSE
      pdots$log <- dots$log
      pdots$adj <- par("adj")

      ####

      co.x <- fct(cutoff.x)(data.x)
      co.y <- fct(cutoff.y)(data.y)
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

      do.call(plot, args = c(list(x = ndata.x,ndata.y, type = "p"), pdots))
      do.call(box,args=c(adots))
      do.call(abline, args = c(list(h=co.y), abdots))
      do.call(abline, args = c(list(v=co.x), abdots))
      if(length(id.xy))
         do.call(text, args = c(list(ndata.x[id.xy], ndata.y[id.xy],
                                labels=lab.pts[id.xy]), tdots))
      return(list(id.x=id0.x, id.y= id0.y, id.xy = id0.xy,
             qtx = quantile(ndata.x), qty = quantile(ndata.y),
             cutoff.x.v = co.x, cutoff.y.v = co.y
             ))
}
