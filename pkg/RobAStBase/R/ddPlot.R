setMethod("ddPlot", signature = signature(data = "matrix"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff){
       mc <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(1)))
       dots <- mc$"..."

do.call(.ddPlot.MatNtNtCoCo, args = c(list(data=data), dots, list(dist.x = mc$dist.x,
       dist.y = mc$dist.y, cutoff.x = mc$cutoff.x, cutoff.y = mc$cutoff.y,
       cutoff.quantile.x = mc$cutoff.quantile.x, cutoff.quantile.y = mc$cutoff.quantile.y,
       transform.x  = mc$transform.x, transform.y = mc$transform.y,
       id.n = mc$id.n, lab.pts = mc$lab.pts, adj = mc$adj, cex.idn = mc$cex.idn,
       col.idn = mc$col.idn, lty.cutoff = mc$lty.cutoff,
       lwd.cutoff = mc$lwd.cutoff, col.cutoff = mc$col.cutoff)))
})

setMethod("ddPlot", signature = signature(data = "data.frame"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff){

         mc <- match.call(call = sys.call(sys.parent(1)))
         mc$data <- t(as.matrix(data))
         do.call(ddPlot,args = as.list(mc[-1]))
       })

setMethod("ddPlot", signature = signature(data = "numeric"),
          function(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff){

         mc <- match.call(call = sys.call(sys.parent(1)))
         mc$data <- matrix(data,nrow=1)
         do.call(ddPlot,args = as.list(mc[-1]))
       })
