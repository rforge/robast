setMethod("ddPlot", signature = signature(data = "matrix"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, cex.pts = 1,lab.pts, jitter.pts = 0, alpha.trsp = NA, adj =0, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff,
       text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jitter.fac, jitter.tol = .Machine$double.eps,doplot = TRUE){

       if(missing(dist.x)) dist.x <- NormType()
       if(missing(dist.y)) dist.y <- NormType()
       if(missing(cutoff.x)) cutoff.x <- NULL
       if(missing(cutoff.y)) cutoff.y <- NULL
       if(missing(transform.x)) transform.x <- NULL
       if(missing(transform.y)) transform.y <- NULL
       if(missing(id.n)) id.n <-  NULL
       if(missing(lab.pts)) lab.pts <- NULL
       if(missing(cex.idn)) cex.idn <- NULL
       if(missing(col.idn)) col.idn <- NULL
       if(missing(lty.cutoff)) lty.cutoff <- NULL
       if(missing(lwd.cutoff)) lwd.cutoff <- NULL
       if(missing(col.cutoff)) col.cutoff <- NULL
       if(missing(col.abline)) col.abline <- NULL
       if(missing(jitter.fac)) jitter.fac <- NULL

       args0 <- list(data = data,
                     dist.x = dist.x,
                     dist.y = dist.y,
                     cutoff.x = cutoff.x,
                     cutoff.y = cutoff.y,
                     cutoff.quantile.x = cutoff.quantile.x,
                     cutoff.quantile.y = cutoff.quantile.y,
                     transform.x = transform.x,
                     transform.y = transform.y,
                     id.n = id.n,
                     cex.pts = cex.pts,
                     lab.pts = lab.pts,
                     jitter.pts = jitter.pts, alpha.trsp = alpha.trsp, adj = adj,
                     cex.idn = cex.idn,
                     col.idn = col.idn,
                     lty.cutoff = lty.cutoff,
                     lwd.cutoff = lwd.cutoff,
                     col.cutoff = col.cutoff,
                     text.abline = text.abline, text.abline.x = text.abline.x,
                     text.abline.y = text.abline.y, cex.abline = cex.abline,
                     col.abline = col.abline,
                     font.abline = font.abline,
                     adj.abline = adj.abline, text.abline.x.x = text.abline.x.x,
                     text.abline.x.y = text.abline.x.y,
                     text.abline.y.x = text.abline.y.x,
                     text.abline.y.y = text.abline.y.y,
                     text.abline.x.fmt.cx = text.abline.x.fmt.cx,
                     text.abline.x.fmt.qx = text.abline.x.fmt.cx,
                     text.abline.y.fmt.cy = text.abline.y.fmt.cy,
                     text.abline.y.fmt.qy = text.abline.y.fmt.qy,
                     jitter.fac = jitter.fac,
                     jitter.tol = jitter.tol,
                     doplot = doplot)

       mc <- match.call(expand.dots = TRUE, call = sys.call(sys.parent(1)))
       dots <- mc$"..."
       plotInfo <- list(call = mc, dots=dots, args=args0)
       mc <- as.list(mc)[-1]
       mc$data <- data
#       ret <- do.call(.ddPlot.MatNtNtCoCo, args = mc)
       ret <- do.call(RobAStBase:::.ddPlot.MatNtNtCoCo, args = mc)
       if(!doplot) return(ret)
       ret$call <- ret$dots <- ret$args <- NULL
       plotInfo <- c(plotInfo, ret)
       class(plotInfo) <- c("plotInfo","DiagnInfo")
       return(invisible(plotInfo))

})

setMethod("ddPlot", signature = signature(data = "data.frame"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, cex.pts = 1,lab.pts, jitter.pts = 0, alpha.trsp = NA, adj =0, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jitter.fac, jitter.tol = .Machine$double.eps,doplot = TRUE){

       if(missing(dist.x)) dist.x <- NormType()
       if(missing(dist.y)) dist.y <- NormType()
       if(missing(cutoff.x)) cutoff.x <- NULL
       if(missing(cutoff.y)) cutoff.y <- NULL
       if(missing(transform.x)) transform.x <- NULL
       if(missing(transform.y)) transform.y <- NULL
       if(missing(id.n)) id.n <-  NULL
       if(missing(lab.pts)) lab.pts <- NULL
       if(missing(cex.idn)) cex.idn <- NULL
       if(missing(col.idn)) col.idn <- NULL
       if(missing(lty.cutoff)) lty.cutoff <- NULL
       if(missing(lwd.cutoff)) lwd.cutoff <- NULL
       if(missing(col.cutoff)) col.cutoff <- NULL
       if(missing(col.abline)) col.abline <- NULL
       if(missing(jitter.fac)) jitter.fac <- NULL

       args0 <- list(data = data,
                     dist.x = dist.x,
                     dist.y = dist.y,
                     cutoff.x = cutoff.x,
                     cutoff.y = cutoff.y,
                     cutoff.quantile.x = cutoff.quantile.x,
                     cutoff.quantile.y = cutoff.quantile.y,
                     transform.x = transform.x,
                     transform.y = transform.y,
                     id.n = id.n,
                     cex.pts = cex.pts,
                     lab.pts = lab.pts,
                     jitter.pts = jitter.pts, alpha.trsp = alpha.trsp, adj = adj,
                     cex.idn = cex.idn,
                     col.idn = col.idn,
                     lty.cutoff = lty.cutoff,
                     lwd.cutoff = lwd.cutoff,
                     col.cutoff = col.cutoff,
                     text.abline = text.abline, text.abline.x = text.abline.x,
                     text.abline.y = text.abline.y, cex.abline = cex.abline,
                     col.abline = col.abline,
                     font.abline = font.abline,
                     adj.abline = adj.abline, text.abline.x.x = text.abline.x.x,
                     text.abline.x.y = text.abline.x.y,
                     text.abline.y.x = text.abline.y.x,
                     text.abline.y.y = text.abline.y.y,
                     text.abline.x.fmt.cx = text.abline.x.fmt.cx,
                     text.abline.x.fmt.qx = text.abline.x.fmt.cx,
                     text.abline.y.fmt.cy = text.abline.y.fmt.cy,
                     text.abline.y.fmt.qy = text.abline.y.fmt.qy,
                     jitter.fac = jitter.fac,
                     jitter.tol = jitter.tol,
                     doplot = doplot)

       mc <- match.call(expand.dots = TRUE, call = sys.call(sys.parent(1)))
       dots <- mc$"..."
       plotInfo <- list(call = mc, dots=dots, args=args0)

       mc$data <- t(as.matrix(data))
       ret <- do.call(ddPlot,args = as.list(mc[-1]))
       if(!doplot) return(ret)
       ret$call <- ret$dots <- ret$args <- NULL
       plotInfo <- c(plotInfo, ret)
       class(plotInfo) <- c("plotInfo","DiagnInfo")
       return(invisible(plotInfo))
       })

setMethod("ddPlot", signature = signature(data = "numeric"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, cex.pts = 1,lab.pts, jitter.pts = 0, alpha.trsp = NA, adj =0, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff,
       text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jitter.fac, jitter.tol = .Machine$double.eps, doplot = TRUE){

       if(missing(dist.x)) dist.x <- NormType()
       if(missing(dist.y)) dist.y <- NormType()
       if(missing(cutoff.x)) cutoff.x <- NULL
       if(missing(cutoff.y)) cutoff.y <- NULL
       if(missing(transform.x)) transform.x <- NULL
       if(missing(transform.y)) transform.y <- NULL
       if(missing(id.n)) id.n <-  NULL
       if(missing(lab.pts)) lab.pts <- NULL
       if(missing(cex.idn)) cex.idn <- NULL
       if(missing(col.idn)) col.idn <- NULL
       if(missing(lty.cutoff)) lty.cutoff <- NULL
       if(missing(lwd.cutoff)) lwd.cutoff <- NULL
       if(missing(col.cutoff)) col.cutoff <- NULL
       if(missing(col.abline)) col.abline <- NULL
       if(missing(jitter.fac)) jitter.fac <- NULL

       args0 <- list(data = data,
                     dist.x = dist.x,
                     dist.y = dist.y,
                     cutoff.x = cutoff.x,
                     cutoff.y = cutoff.y,
                     cutoff.quantile.x = cutoff.quantile.x,
                     cutoff.quantile.y = cutoff.quantile.y,
                     transform.x = transform.x,
                     transform.y = transform.y,
                     id.n = id.n,
                     cex.pts = cex.pts,
                     lab.pts = lab.pts,
                     jitter.pts = jitter.pts, alpha.trsp = alpha.trsp, adj = adj,
                     cex.idn = cex.idn,
                     col.idn = col.idn,
                     lty.cutoff = lty.cutoff,
                     lwd.cutoff = lwd.cutoff,
                     col.cutoff = col.cutoff,
                     text.abline = text.abline, text.abline.x = text.abline.x,
                     text.abline.y = text.abline.y, cex.abline = cex.abline,
                     col.abline = col.abline,
                     font.abline = font.abline,
                     adj.abline = adj.abline, text.abline.x.x = text.abline.x.x,
                     text.abline.x.y = text.abline.x.y,
                     text.abline.y.x = text.abline.y.x,
                     text.abline.y.y = text.abline.y.y,
                     text.abline.x.fmt.cx = text.abline.x.fmt.cx,
                     text.abline.x.fmt.qx = text.abline.x.fmt.cx,
                     text.abline.y.fmt.cy = text.abline.y.fmt.cy,
                     text.abline.y.fmt.qy = text.abline.y.fmt.qy,
                     jitter.fac = jitter.fac,
                     jitter.tol = jitter.tol,
                     doplot = doplot)

       mc <- match.call(expand.dots = TRUE, call = sys.call(sys.parent(1)))
       dots <- mc$"..."
       plotInfo <- list(call = mc, dots=dots, args=args0)

       mc$data <- matrix(data,nrow=1)
       ret <- do.call(ddPlot,args = as.list(mc[-1]))
       if(!doplot) return(ret)
       ret$call <- ret$dots <- ret$args <- NULL
       plotInfo <- c(plotInfo, ret)
       class(plotInfo) <- c("plotInfo","DiagnInfo")
       return(invisible(plotInfo))
       })
