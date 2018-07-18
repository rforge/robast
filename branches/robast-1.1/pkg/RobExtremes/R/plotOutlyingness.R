##########################################
##                                      ## 
##    Wrapper for outlyingnessPlot.R    ##
##                                      ##
##                                      ##
##########################################

if(FALSE){

##@x - dataset
##@X - random variable
##@fam - parameter family
##@alpha - confidence level for quantile
#
plotOutlyingness = function(x,alpha=0.99,X=GPareto(),fam=GParetoFamily(),
                           dist.y =  NormType(),
                           cutoff.x = cutoff(), cutoff.y = cutoff.sememp(),
                           ...,
                           id.n,
                           cex.pts = 1,
                           lab.pts,
                           jitter.pts = 0,
                           alpha.trsp = NA,
                           adj = 0.1,
                           pch = 16,
                           cex = 1.5,
                           col = rgb(152,152,152,maxColorValue=255),
                           cex.idn = 1.7,
                           col.idn = rgb(102,102,102,maxColorValue=255),
                           lty.cutoff,
                           lwd.cutoff = 3,
                           col.cutoff = rgb(202,202,202,maxColorValue=255),
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
                           robCov.x = TRUE,
                           robCov.y = TRUE,
                           tf.x = function(x)apply(x,2,log),
                           tf.y = function(x)x[1,],
                           jitter.fac = 10,
                           jitter.tol=.Machine$double.eps,
                           cex.lab = 1.5,
                           col.lab="red",
                           main = "Outlyingness Plot for Extreme Value Distributions",
                           cex.main = 1.5,
                           xlab="Theoretical log-quantiles",
                           ylab="Mahalanobis distance",
                           doplot = TRUE
                           ){


  mc <- match.call()
  dots <- match.call(expand.dots = FALSE)$"..."
  args0 <- list(x = x, alpha = alpha, X = X, fam = fam, dist.y = dist.y,
               cutoff.x = cutoff.x, cutoff.y = cutoff.y,
               id.n = if(missing(id.n)) NULL else id.n,
               cex.pts = cex.pts, lab.pts = if(missing(lab.pts)) NULL else lab.pts,
               jitter.pts = jitter.pts, alpha.trsp = alpha.trsp, adj = adj,
               pch = pch, cex = cex, col = col, cex.idn = cex.idn, col.idn = col.idn,
               lty.cutoff = if(missing(lty.cutoff)) NULL else lty.cutoff,
               lwd.cutoff = lwd.cutoff, col.cutoff = col.cutoff,
               text.abline =  text.abline, text.abline.x = text.abline.x,
               text.abline.y = text.abline.y, cex.abline = cex.abline,
               col.abline = if(missing(col.abline)) col.cutoff else col.abline,
               font.abline = font.abline, adj.abline = adj.abline,
               text.abline.x.x = text.abline.x.x, text.abline.x.y = text.abline.x.y,
               text.abline.y.x = text.abline.y.x, text.abline.y.y = text.abline.y.y,
               text.abline.x.fmt.cx = text.abline.x.fmt.cx,
               text.abline.x.fmt.qx = text.abline.x.fmt.qx,
               text.abline.y.fmt.cy = text.abline.y.fmt.cy,
               text.abline.y.fmt.qy = text.abline.y.fmt.qy, robCov.x = robCov.x,
               robCov.y = robCov.y, tf.x = tf.x, tf.y = tf.y,
               jitter.fac = jitter.fac, jitter.tol = jitter.tol, cex.lab = cex.lab,
               col.lab = col.lab, main = main, cex.main = cex.main, xlab = xlab,
               ylab = ylab, doplot = doplot)
  plotInfo <- list(call = mc, dots=dots, args=args0)


  ##projection distance
  qfun = function(x){p0 = p(X)(x); q0 = q.l(X)(p0)}
  QProj <- function(){new("NormType", name="Quantiles", fct=qfun)}

  ##logarithmic representation (for distributions with positive support)
  fam@distribution = log(fam@distribution)

  ##classical IC
  ICmle <- optIC(model=fam,risk=asCov())

  ##parameter for plotting
  par(cex=1,bty="n")

  ##call of routine from RobAStBase
  plotInfo$outlyingPlotICArgs <- c(list(data = x, IC.x = ICmle, IC.y = ICmle,
      dist.x = QProj(), dist.y = dist.y, cutoff.x = cutoff.x, cutoff.y = cutoff.y),
      dots,list(cutoff.quantile.y = cutoff.quantile.x,
        cutoff.quantile.x = cutoff.quantile.y, id.n = id.n, cex.pts = cex.pts,
        lab.pts = lab.pts, jitter.pts = jitter.pts, alpha.trsp = alpha.trsp,
        adj = adj, pch = pch, cex = cex, col = col, cex.idn = cex.idn,
        col.idn = col.idn, lty.cutoff = lty.cutoff, lwd.cutoff = lwd.cutoff,
        col.cutoff = col.cutoff, text.abline =  text.abline,
        text.abline.x = text.abline.x, text.abline.y = text.abline.y,
        cex.abline = cex.abline, col.abline = col.abline,
        font.abline = font.abline, adj.abline = adj.abline,
        text.abline.x.x = text.abline.x.x, text.abline.x.y = text.abline.x.y,
        text.abline.y.x = text.abline.y.x, text.abline.y.y = text.abline.y.y,
        text.abline.x.fmt.cx = text.abline.x.fmt.cx,
        text.abline.x.fmt.qx = text.abline.x.fmt.qx,
        text.abline.y.fmt.cy = text.abline.y.fmt.cy,
        text.abline.y.fmt.qy = text.abline.y.fmt.qy, robCov.x = robCov.x,
        robCov.y = robCov.y, tf.x = tf.x, tf.y = tf.y, jitter.fac = jitter.fac,
        jitter.tol = jitter.tol, cex.lab = cex.lab, col.lab = col.lab,
        main = main, cex.main = cex.main, xlab = xlab, ylab = ylab, doplot = doplot))
retV <- do.call(outlyingPlotIC,args=plotInfo$outlyingPlotICArgs)
retV$args <- NULL
retV$dots <- NULL
retV$call <- NULL
plotInfo <- c(PlotInfo, retV)
class(plotInfo) <- c("plotInfo","DiagnInfo")
return(invisible(plotInfo))
}

##Example
X= GPareto()
fam = GParetoFamily()
x = r(X)(1000)
plotOutlyingness(x,alpha=0.9,X=X,fam=fam)
}
