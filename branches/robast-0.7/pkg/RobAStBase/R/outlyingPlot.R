outlyingPlotIC <- function(data, IC.x, IC.y = IC.x, dist.x = NormType(),
                         dist.y, cutoff.y = cutoff.chisq(), cutoff.x = cutoff.sememp(), ...,
                         cutoff.quantile.x = 0.95,
                         cutoff.quantile.y = cutoff.quantile.x,
                         id.n, lab.pts, adj, cex.idn,
                         col.idn, lty.cutoff, lwd.cutoff, col.cutoff,
                         main = gettext("Outlyingness by means of a distance-distance plot")
                         ){
     mc <- as.list(match.call(expand.dots = FALSE))[-1]
     dots <- mc$"..."
     if(is.null(dots$xlim)) dots$xlim <- TRUE
     if(is.null(dots$ylim)) dots$ylim <- TRUE
     if(is.null(mc$cutoff.quantile.x)) mc$cutoff.quantile.x <- 0.95
     if(is.null(mc$cutoff.quantile.y)) mc$cutoff.quantile.y <- cutoff.quantile.x
     if(is.null(mc$cutoff.x)) mc$cutoff.x <- cutoff.sememp()
     if(is.null(mc$cutoff.y)) mc$cutoff.y <- cutoff.chisq()
     if(missing(IC.x)) stop("Argument 'IC.x' must be given as argument to 'outlyingPlot'")
     if(missing(data)) stop("Argument 'data' must be given as argument to 'outlyingPlot'")

     if(missing(dist.y)){
        if("asCov" %in% names(Risks(IC.y)))
            if(is.matrix(Risks(IC.y)$asCov) || length(Risks(IC.y)$asCov) == 1)
               asVar <- Risks(IC.y)$asCov
            else
               asVar <- Risks(IC.y)$asCov$value
        else
            asVar <- getRiskIC(IC.y, risk = asCov())$asCov$value

        asVar <- PosSemDefSymmMatrix(solve(asVar))
        mc$dist.y <- QFNorm(name = gettext("Mahalonobis-Norm"), QuadForm = asVar)
     }
     if(missing(dist.x))
        mc$dist.x <- NormType()

     tf.x <- function(x) apply(x,2,function(xx) evalIC(IC.x,xx))
     tf.y <- function(x) apply(x,2,function(xx) evalIC(IC.y,xx))


     do.call(ddPlot,args=c(list(data=data), dots, list(dist.x = mc$dist.x,
       dist.y = mc$dist.y, cutoff.x = mc$cutoff.x, cutoff.y = mc$cutoff.y,
       cutoff.quantile.x = mc$cutoff.quantile.x, cutoff.quantile.y = mc$cutoff.quantile.y,
       transform.x = tf.x, transform.y = tf.y,
       id.n = mc$id.n, lab.pts = mc$lab.pts, adj = mc$adj, cex.idn = mc$cex.idn,
       col.idn = mc$col.idn, lty.cutoff = mc$lty.cutoff,
       lwd.cutoff = mc$lwd.cutoff, col.cutoff = mc$col.cutoff, main = main)))

     }

