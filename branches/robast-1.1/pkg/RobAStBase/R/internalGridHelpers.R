.getDimsTD <- function(L2Fam,to.draw.arg){
  trafO <- trafo(L2Fam@param)
  dims  <- nrow(trafO)
  dimnms  <- c(rownames(trafO))
  if(is.null(dimnms))
     dimnms <- paste("dim",1:dims,sep="")
  to.draw <- 1:dims
  if(! is.null(to.draw.arg)){
       if(is.character(to.draw.arg))
          to.draw <- pmatch(to.draw.arg, dimnms)
       else if(is.numeric(to.draw.arg))
               to.draw <- to.draw.arg
  }
  return(length(to.draw))
}


.producePanelFirstS <- function(panelFirst,IC,to.draw.arg, isInfoPlot=FALSE,
                                x.ticks, scaleX, scaleX.fct,
                                y.ticks, scaleY, scaleY.fct){


  L2Fam <- eval(IC@CallL2Fam)
  if(is.null(scaleX.fct)) scaleX.fct <- p(L2Fam)
  ndim <- .getDimsTD(L2Fam,to.draw.arg)
  if(is.null(scaleY.fct)){
     scaleY.fct <- .fillList(pnorm,ndim)
  }else{
     scaleY.fct <- .fillList(scaleY.fct,ndim)
  }
  ..y.ticks <- .fillList(y.ticks,ndim)
  .xticksS <- substitute({
            .x.ticks <- x.ticks0
            if(is.null(.x.ticks))
               .x.ticks <- axTicks(1, axp=par("xaxp"), usr=par("usr"))
            scaleX00 <- FALSE
            if(!is.null(scaleX0)) scaleX00 <- scaleX0
            if(scaleX00) .x.ticks <- scaleX.fct0(.x.ticks)
            },
            list(x.ticks0 = x.ticks, scaleX0 = scaleX, scaleX.fct0 = scaleX.fct)
            )

  getYI <- if(isInfoPlot){
      substitute({
             i0 <- if(exists("i")) get("i") else 1
            .y.ticks <- if(.absInd) NULL else .y.ticksL[[i0]]
      })
  }else{
      substitute({
            .y.ticks <- .y.ticksL[[i]]
      })

  }

  assYI <- if(isInfoPlot){
      substitute({
             i0 <- if(exists("i")) get("i") else 1
             if(.absInd) .y.ticks <- scaleY.fct0[[i0]](.y.ticks)
      }, list(scaleY.fct0 = scaleY.fct))
  }else{
      substitute({
            .y.ticks <- scaleY.fct0[[i]](.y.ticks)
      }, list(scaleY.fct0 = scaleY.fct))

  }

  .yticksS <- substitute({
            .y.ticksL <- y.ticks0
            getYI0
            if(is.null(.y.ticks))
               .y.ticks <- axTicks(2, axp=par("yaxp"), usr=par("usr"))
            scaleY00 <- FALSE
            if(!is.null(scaleY0)) scaleY00 <- scaleY0
            if(scaleY00) assYI0
            },
            list(y.ticks0 = y.ticks, scaleY0 = scaleY, scaleY.fct0 = scaleY.fct,
                 getYI0 = getYI, assYI0 = assYI)
            )

  ..panelFirst  <- panelFirst
  if(length(panelFirst)){
     for(i in 1:length(panelFirst)){
         ..panelFirst[[i]] <- substitute({
             pFi
             .xticksS0
             .yticksS0
             abline(v=.x.ticks,col= "lightgray",
                    lty = "dotted", lwd = par("lwd"))
             abline(h=.y.ticks,col= "lightgray",
                    lty = "dotted", lwd = par("lwd"))
          },list(pFi = if(is.null(panelFirst[[i]])) expression({}) else panelFirst[[i]],
                 .xticksS0 = .xticksS, .yticksS0 = .yticksS)
          )
     }
   }
   return(..panelFirst)
}
