setMethod("ddPlot", signature = signature(data = "matrix"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, cex.pts = 1,lab.pts, jit.pts = 0, alpha.trsp = NA, adj =0, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jit.fac, jit.tol = .Machine$double.eps,doplot = TRUE){

       args0 <- list(data = data,
                     dist.x = if(!missing(dist.x)) dist.x else NULL,
                     dist.y = if(!missing(dist.y)) dist.y else NULL,
                     cutoff.x = if(!missing(cutoff.x)) cutoff.x else NULL,
                     cutoff.y = if(!missing(cutoff.y)) cutoff.y else NULL,
                     cutoff.quantile.x = cutoff.quantile.x,
                     cutoff.quantile.y = cutoff.quantile.y,
                     transform.x = if(!missing(transform.x)) transform.x else NULL,
                     transform.y = if(!missing(transform.y)) transform.y else NULL,
                     id.n = if(!missing(id.n)) id.n else NULL,
                     cex.pts = cex.pts,
                     lab.pts = if(!missing(lab.pts)) lab.pts else NULL,
                     jit.pts = jit.pts, alpha.trsp = alpha.trsp, adj = adj,
                     cex.idn =if(!missing(cex.idn)) cex.idn else NULL,
                     col.idn =if(!missing(col.idn)) col.idn else NULL,
                     lty.cutoff =if(!missing(lty.cutoff)) lty.cutoff else NULL,
                     lwd.cutoff =if(!missing(lwd.cutoff)) lwd.cutoff else NULL,
                     col.cutoff =if(!missing(col.cutoff)) col.cutoff else NULL,
                     text.abline = text.abline, text.abline.x = text.abline.x,
                     text.abline.y = text.abline.y, cex.abline = cex.abline,
                     col.abline = if(!missing(col.abline)) col.abline else NULL,
                     font.abline = font.abline,
                     adj.abline = adj.abline, text.abline.x.x = text.abline.x.x,
                     text.abline.x.y = text.abline.x.y,
                     text.abline.y.x = text.abline.y.x,
                     text.abline.y.y = text.abline.y.y,
                     text.abline.x.fmt.cx = text.abline.x.fmt.cx,
                     text.abline.x.fmt.qx = text.abline.x.fmt.cx,
                     text.abline.y.fmt.cy = text.abline.y.fmt.cy,
                     text.abline.y.fmt.qy = text.abline.y.fmt.qy,
                     jit.fac = if(!missing(jit.fac)) jit.fac else NULL,
                     jit.tol = jit.tol,
                     doplot = doplot)

       mc <- match.call(expand.dots = TRUE, call = sys.call(sys.parent(1)))
       dots <- mc$"..."
       plotInfo <- list(call = mc, dots=dots, args=args0)
       mc <- as.list(mc)[-1]
       mc$data <- data
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
       id.n, cex.pts = 1,lab.pts, jit.pts = 0, alpha.trsp = NA, adj =0, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jit.fac, jit.tol = .Machine$double.eps,doplot = TRUE){

       args0 <- list(data = data,
                     dist.x = if(!missing(dist.x)) dist.x else NULL,
                     dist.y = if(!missing(dist.y)) dist.y else NULL,
                     cutoff.x = if(!missing(cutoff.x)) cutoff.x else NULL,
                     cutoff.y = if(!missing(cutoff.y)) cutoff.y else NULL,
                     cutoff.quantile.x = cutoff.quantile.x,
                     cutoff.quantile.y = cutoff.quantile.y,
                     transform.x = if(!missing(transform.x)) transform.x else NULL,
                     transform.y = if(!missing(transform.y)) transform.y else NULL,
                     id.n = if(!missing(id.n)) id.n else NULL,
                     cex.pts = cex.pts,
                     lab.pts = if(!missing(lab.pts)) lab.pts else NULL,
                     jit.pts = jit.pts, alpha.trsp = alpha.trsp, adj = adj,
                     cex.idn =if(!missing(cex.idn)) cex.idn else NULL,
                     col.idn =if(!missing(col.idn)) col.idn else NULL,
                     lty.cutoff =if(!missing(lty.cutoff)) lty.cutoff else NULL,
                     lwd.cutoff =if(!missing(lwd.cutoff)) lwd.cutoff else NULL,
                     col.cutoff =if(!missing(col.cutoff)) col.cutoff else NULL,
                     text.abline = text.abline, text.abline.x = text.abline.x,
                     text.abline.y = text.abline.y, cex.abline = cex.abline,
                     col.abline = if(!missing(col.abline)) col.abline else NULL,
                     font.abline = font.abline,
                     adj.abline = adj.abline, text.abline.x.x = text.abline.x.x,
                     text.abline.x.y = text.abline.x.y,
                     text.abline.y.x = text.abline.y.x,
                     text.abline.y.y = text.abline.y.y,
                     text.abline.x.fmt.cx = text.abline.x.fmt.cx,
                     text.abline.x.fmt.qx = text.abline.x.fmt.cx,
                     text.abline.y.fmt.cy = text.abline.y.fmt.cy,
                     text.abline.y.fmt.qy = text.abline.y.fmt.qy,
                     jit.fac = if(!missing(jit.fac)) jit.fac else NULL,
                     jit.tol = jit.tol,
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
       id.n, cex.pts = 1,lab.pts, jit.pts = 0, alpha.trsp = NA, adj =0, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff,
       text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL,
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "%7.2f", text.abline.x.fmt.qx = "%4.2f%%",
       text.abline.y.fmt.cy = "%7.2f", text.abline.y.fmt.qy = "%4.2f%%",
       jit.fac, jit.tol = .Machine$double.eps, doplot = TRUE){

       args0 <- list(data = data,
                     dist.x = if(!missing(dist.x)) dist.x else NULL,
                     dist.y = if(!missing(dist.y)) dist.y else NULL,
                     cutoff.x = if(!missing(cutoff.x)) cutoff.x else NULL,
                     cutoff.y = if(!missing(cutoff.y)) cutoff.y else NULL,
                     cutoff.quantile.x = cutoff.quantile.x,
                     cutoff.quantile.y = cutoff.quantile.y,
                     transform.x = if(!missing(transform.x)) transform.x else NULL,
                     transform.y = if(!missing(transform.y)) transform.y else NULL,
                     id.n = if(!missing(id.n)) id.n else NULL,
                     cex.pts = cex.pts,
                     lab.pts = if(!missing(lab.pts)) lab.pts else NULL,
                     jit.pts = jit.pts, alpha.trsp = alpha.trsp, adj = adj,
                     cex.idn =if(!missing(cex.idn)) cex.idn else NULL,
                     col.idn =if(!missing(col.idn)) col.idn else NULL,
                     lty.cutoff =if(!missing(lty.cutoff)) lty.cutoff else NULL,
                     lwd.cutoff =if(!missing(lwd.cutoff)) lwd.cutoff else NULL,
                     col.cutoff =if(!missing(col.cutoff)) col.cutoff else NULL,
                     text.abline = text.abline, text.abline.x = text.abline.x,
                     text.abline.y = text.abline.y, cex.abline = cex.abline,
                     col.abline = if(!missing(col.abline)) col.abline else NULL,
                     font.abline = font.abline,
                     adj.abline = adj.abline, text.abline.x.x = text.abline.x.x,
                     text.abline.x.y = text.abline.x.y,
                     text.abline.y.x = text.abline.y.x,
                     text.abline.y.y = text.abline.y.y,
                     text.abline.x.fmt.cx = text.abline.x.fmt.cx,
                     text.abline.x.fmt.qx = text.abline.x.fmt.cx,
                     text.abline.y.fmt.cy = text.abline.y.fmt.cy,
                     text.abline.y.fmt.qy = text.abline.y.fmt.qy,
                     jit.fac = if(!missing(jit.fac)) jit.fac else NULL,
                     jit.tol = jit.tol,
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
