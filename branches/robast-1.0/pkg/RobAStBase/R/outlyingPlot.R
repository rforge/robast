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
                           jitt.pts = 0,
                           alpha.trsp = NA,
                           adj,
                           cex.idn,
                           col.idn, 
                           lty.cutoff, 
                           lwd.cutoff, 
                           col.cutoff,
                           robCov.x = TRUE,
                           robCov.y = TRUE,
                           tf.x = data,
                           tf.y = data,
                           jitt.fac=10,
                           doplot = TRUE,
                           main = gettext("Outlyingness \n by means of a distance-distance plot")
                           ){
     mc <- as.list(match.call(expand.dots = FALSE))[-1]
     dots <- mc$"..."

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
         CMcd <- PosSemDefSymmMatrix(getCov(CovMcd(devIC,alpha=0.5)))
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
            CMcd <- PosSemDefSymmMatrix(getCov(CovMcd(devIC,alpha=0.5)))
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


    if(missing(tf.x)){
     tf.x <- function(x) apply(x,2,function(xx) evalIC(IC.x,xx))
     }else{tf.x <- mc$tf.x}
    if(missing(tf.y)){
     tf.y <- function(x) apply(x,2,function(xx) evalIC(IC.y,xx))
     }else{tf.y <- mc$tf.y}

    if(!missing(cutoff.x)) assign("..ICloc", IC.x, envir=environment(fct(cutoff.x)))
    if(!missing(cutoff.y)) assign("..ICloc", IC.y, envir=environment(fct(cutoff.y)))
     do.call(ddPlot,args=c(list(data=data),dots,
       list(dist.x = mc$dist.x,
       dist.y = mc$dist.y, 
       cutoff.x = cutoff.x,
       cutoff.y = cutoff.y,
       cutoff.quantile.x = mc$cutoff.quantile.x, 
       cutoff.quantile.y = mc$cutoff.quantile.y,
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
       jitt.fac = mc$jitt.fac,
       doplot = doplot,
       main = main)))

     }

