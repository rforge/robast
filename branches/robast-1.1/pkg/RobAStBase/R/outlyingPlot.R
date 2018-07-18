outlyingPlotIC <- function(data, 
                           IC.x, 
                           IC.y = IC.x,
                           dist.x = NormType(),
                           dist.y, 
                           cutoff.x = cutoff.sememp(0.95),
                           cutoff.y = cutoff.chisq(0.95), 
                           ...,
                           cutoff.quantile.x = 0.95,
                           cutoff.quantile.y = cutoff.quantile.x,
                           id.n,
                           cex.pts = 1,
                           lab.pts,
                           jitter.pts = 0,
                           alpha.trsp = NA,
                           adj,
                           cex.idn,
                           col.idn, 
                           lty.cutoff, 
                           lwd.cutoff, 
                           col.cutoff,
                           text.abline = TRUE,
                           text.abline.x = NULL,
                           text.abline.y = NULL,
                           cex.abline = par("cex"),
                           col.abline = col.cutoff,
                           font.abline = par("font"),
                           adj.abline = c(0,0),
                           text.abline.x.x = NULL,
                           text.abline.x.y = NULL,
                           text.abline.y.x = NULL,
                           text.abline.y.y = NULL,
                           text.abline.x.fmt.cx = "%7.2f",
                           text.abline.x.fmt.qx = "%4.2f%%",
                           text.abline.y.fmt.cy = "%7.2f",
                           text.abline.y.fmt.qy = "%4.2f%%",
                           robCov.x = TRUE,
                           robCov.y = TRUE,
                           tf.x = NULL,
                           tf.y = NULL,
                           jitter.fac=10,
                           jitter.tol=.Machine$double.eps,
                           doplot = TRUE,
                           main = gettext("Outlyingness \n by means of a distance-distance plot")
                           ){

        if(missing(dist.x)) dist.x <- NormType()
        if(missing(dist.y)) dist.y <- NULL
        if(missing(id.n)) id.n <- NULL
        if(missing(lab.pts)) lab.pts <- NULL
        if(missing(adj)) adj <- NULL
        if(missing(cex.idn)) cex.idn <- NULL
        if(missing(col.idn)) col.idn <- NULL
        if(missing(lty.cutoff)) lty.cutoff <- NULL
        if(missing(lwd.cutoff)) lwd.cutoff <- NULL
        if(missing(col.cutoff)) col.cutoff <- NULL

        args0 <- list(data = data, IC.x = IC.x, IC.y = IC.y,
                      dist.x = dist.x, dist.y = dist.y,
                      cutoff.x = cutoff.x, cutoff.y = cutoff.y,
                      cutoff.quantile.x = cutoff.quantile.x,
                      cutoff.quantile.y = cutoff.quantile.y,
                      id.n = id.n, cex.pts = cex.pts, lab.pts = lab.pts,
                      jitter.pts = jitter.pts, alpha.trsp = alpha.trsp,
                      adj = adj, cex.idn = cex.idn, col.idn = col.idn,
                      lty.cutoff = lty.cutoff, lwd.cutoff = lwd.cutoff,
                      col.cutoff = col.cutoff,
                      text.abline =  text.abline,
                      text.abline.x = text.abline.x,
                      text.abline.y = text.abline.y,
                      cex.abline = cex.abline,
                      col.abline = col.abline,
                      font.abline = font.abline,
                      adj.abline = adj.abline,
                      text.abline.x.x = text.abline.x.x,
                      text.abline.x.y = text.abline.x.y,
                      text.abline.y.x = text.abline.y.x,
                      text.abline.y.y = text.abline.y.y,
                      text.abline.x.fmt.cx = text.abline.x.fmt.cx,
                      text.abline.x.fmt.qx = text.abline.x.fmt.qx,
                      text.abline.y.fmt.cy = text.abline.y.fmt.cy,
                      text.abline.y.fmt.qy = text.abline.y.fmt.qy,
                      robCov.x = robCov.x,robCov.y = robCov.x,
                      tf.x = tf.x, tf.y = tf.y, jitter.fac=jitter.fac,
                      jitter.tol = jitter.tol, doplot = doplot,
                      main = main)
     mc <- match.call(expand.dots = FALSE)
     dots <- mc$"..."
     plotInfo <- list(call = mc, dots=dots, args=args0)

     mc1 <- mc[-1]

     if(is.null(dots$xlim)) dots$xlim <- TRUE
     if(is.null(dots$ylim)) dots$ylim <- TRUE
     if(is.null(mc$cutoff.quantile.x)) mc$cutoff.quantile.x <- 0.95
     if(is.null(mc$cutoff.quantile.y)) mc$cutoff.quantile.y <- cutoff.quantile.x
     if(is.null(mc$cutoff.x)) mc$cutoff.x <- cutoff.sememp(mc$cutoff.quantile.x)
     if(is.null(mc$cutoff.y)) mc$cutoff.y <- cutoff.chisq(mc$cutoff.quantile.y)
     if(missing(IC.x)) stop("Argument 'IC.x' must be given as argument to 'outlyingPlot'")
     if(missing(data)) stop("Argument 'data' must be given as argument to 'outlyingPlot'")
  
  if(missing(dist.x)){
        #mc$dist.x <- NormType()
    if(robCov.x){
      evIC = evalIC(IC.x,as.matrix(data))
      if(is.null(dim(evIC))){
         asVar <- PosSemDefSymmMatrix(mad(evIC)^2)
         if(asVar < 1e-8) asVar <- 1
      }else{
         dimevIC <- dim(evIC)[1]
         devIC <- data.frame(t(evIC[1:dimevIC,,drop=FALSE]))
         CMcd <- PosSemDefSymmMatrix(rrcov::getCov(rrcov::CovMcd(devIC,alpha=0.5)))
         asVar <- CMcd
#         asVar <- solve(CMcd)
#         cat("\n", sep="", gettext("Robust asVar"), ":\n")
#         print(asVar)
      }
      #cat("\nRobust asVar:") ;print("KKKKK")
      #print(asVar)
   }else{if("asCov" %in% names(Risks(IC.x)))
      if(is.matrix(Risks(IC.x)$asCov) || length(Risks(IC.x)$asCov) == 1)
               {asVar <- Risks(IC.x)$asCov
                  }
         else
               {asVar <- Risks(IC.x)$asCov$value 
                  }
         else
            {asVar <- getRiskIC(IC.x, risk = asCov())$asCov$value
                   }
       }
    
#       asVar <- PosSemDefSymmMatrix(solve(asVar))
       mc$dist.x <- QFNorm(name = gettext("Mahalonobis-Norm"), QuadForm = PosSemDefSymmMatrix(solve(asVar)))
      }

     if(missing(dist.y)){
       if(robCov.y){
          evIC <- evalIC(IC.y,as.matrix(data))
          if(is.null(dim(evIC))){
             asVar <- PosSemDefSymmMatrix(mad(evIC)^2)
             if(asVar < 1e-8) asVar <- 1
          }else{
            dimevIC <- dim(evIC)[1]
            devIC <- data.frame(t(evIC[1:dimevIC,,drop=FALSE]))
            CMcd <- PosSemDefSymmMatrix(rrcov::getCov(rrcov::CovMcd(devIC,alpha=0.5)))
            asVar <- CMcd
            cat("Fall 1\n\n")
            print(asVar)
          }
       }else{
            if("asCov" %in% names(Risks(IC.y)))
               if(is.matrix(Risks(IC.y)$asCov) || length(Risks(IC.y)$asCov) == 1)
                  {asVar <- Risks(IC.y)$asCov
                   }
                else{asVar <- Risks(IC.y)$asCov$value
                   }
                else{asVar <- getRiskIC(IC.y, risk = asCov())$asCov$value
                   }
            cat("Fall 2\n\n")
            print(asVar)
       }
     
          mc$dist.y <- QFNorm(name = gettext("Mahalonobis-Norm"), 
                              QuadForm =  PosSemDefSymmMatrix(solve(asVar)))
     }


    if(missing(tf.x)||is.null(tf.x)){
     tf.x <- function(x) apply(x,2,function(xx) evalIC(IC.x,xx))
     }else{tf.x <- mc$tf.x}
    if(missing(tf.y)||is.null(tf.y)){
     tf.y <- function(x) apply(x,2,function(xx) evalIC(IC.y,xx))
     }else{tf.y <- mc$tf.y}

    if(!missing(cutoff.x)) assign("..ICloc", IC.x, envir=environment(fct(cutoff.x)))
    if(!missing(cutoff.y)) assign("..ICloc", IC.y, envir=environment(fct(cutoff.y)))
    plotInfo$ddPlotArgs <- c(list(data=data),dots,
       list(dist.x = mc$dist.x,
       dist.y = mc$dist.y,
       cutoff.x = cutoff.x,
       cutoff.y = cutoff.y,
       cutoff.quantile.x = mc$cutoff.quantile.x,
       cutoff.quantile.y = mc$cutoff.quantile.y,
       jitter.pts = mc$jitter.pts,
       transform.x = tf.x,
       transform.y = tf.y,
       id.n = mc$id.n,
       lab.pts = mc$lab.pts,
       alpha.trsp = alpha.trsp,
       cex.pts = cex.pts,
       adj = mc$adj,
       cex.idn = mc$cex.idn,
       col.idn = mc$col.idn,
       lty.cutoff = mc$lty.cutoff,
       lwd.cutoff = mc$lwd.cutoff,
       col.cutoff = mc$col.cutoff,
       text.abline =  mc$text.abline,
       text.abline.x = mc$text.abline.x,
       text.abline.y = mc$text.abline.y,
       cex.abline = mc$cex.abline,
       col.abline = mc$col.abline,
       font.abline = mc$font.abline,
       adj.abline = mc$adj.abline,
       text.abline.x.x = mc$text.abline.x.x,
       text.abline.x.y = mc$text.abline.x.y,
       text.abline.y.x = mc$text.abline.y.x,
       text.abline.y.y = mc$text.abline.y.y,
       text.abline.x.fmt.cx = mc$text.abline.x.fmt.cx,
       text.abline.x.fmt.qx = mc$text.abline.x.fmt.qx,
       text.abline.y.fmt.cy = mc$text.abline.y.fmt.cy,
       text.abline.y.fmt.qy = mc$text.abline.y.fmt.qy,
       jitter.fac = mc$jitter.fac,
       jitter.tol = mc$jitter.tol,
       doplot = doplot,
       main = main))
     ret <- do.call(ddPlot,args=plotInfo$ddPlotArgs)
     if(!doplot) return(ret)
     ret$args<- NULL
     ret$call<- NULL
     ret$dots<- NULL
     plotInfo <- c(plotInfo,ret)
     class(plotInfo) <- c("plotInfo","DiagnInfo")
     return(invisible(plotInfo))
}

