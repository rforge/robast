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
  }else{if(!is.list(scaleY.fct))
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
  if(is.list(panelFirst) && length(panelFirst)>0){
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
   }else{
    if(!is.list(panelFirst)){
         ..panelFirst <- substitute({
             pFi
             .xticksS0
             .yticksS0
             abline(v=.x.ticks,col= "lightgray",
                    lty = "dotted", lwd = par("lwd"))
             abline(h=.y.ticks,col= "lightgray",
                    lty = "dotted", lwd = par("lwd"))
          },list(pFi = if(is.null(panelFirst)) expression({}) else panelFirst,
                 .xticksS0 = .xticksS, .yticksS0 = .yticksS)
          )
     }
   }
   return(..panelFirst)
}


.producePanelFirstSn <- function(panelFirst,
                                x.ticks, scaleX, scaleX.fct,
                                y.ticks, scaleY, scaleY.fct){


  if(is.null(scaleX.fct)) scaleX.fct <- pnorm
  if(is.null(scaleY.fct)) scaleY.fct <- pnorm

  if(!is.null(x.ticks)){
     .xticksS <- substitute({.x.ticks <- x0}, list(x0 = if(scaleX) x.ticks else scaleX.fct(x.ticks)))
  }else{
     if(!scaleX){
        .xticksS <- substitute({
               .x.ticks <- axTicks(1, axp=par("xaxp"), usr=par("usr"))
               })
     }else{
        .xticksS <- substitute({
               .x.ticks <- fct(axTicks(1, axp=par("xaxp"), usr=par("usr")))
                },list(fct=scaleX.fct))
    }
  }

  if(!is.null(y.ticks)){
     .yticksS <- substitute({.y.ticks <- y0}, list(y0 = if(scaleY) y.ticks else scaleY.fct(y.ticks)))
  }else{
     if(!scaleY){
        .yticksS <- substitute({
               .y.ticks <- axTicks(2, axp=par("yaxp"), usr=par("usr"))
               })
     }else{
        .yticksS <- substitute({
               .y.ticks <- fct(axTicks(2, axp=par("yaxp"), usr=par("usr")))
                },list(fct=scaleY.fct))
     }
  }
  ..panelFirst <- substitute({
        pFi
        .xticksS0
        .yticksS0
        abline(v=.x.ticks,col= "lightgray",
                    lty = "dotted", lwd = par("lwd"))
        abline(h=.y.ticks,col= "lightgray",
                    lty = "dotted", lwd = par("lwd"))
         },list(pFi = if(is.null(panelFirst)) expression({}) else panelFirst,
                 .xticksS0 = .xticksS, .yticksS0 = .yticksS))
  return(..panelFirst)
}

.getX.vec <- function(distr, dims0, lty, x.vec, scaleX, scaleX.fct, scaleX.inv, xm, xM){
        if(!is.null(x.vec)){
           x.vec <- .fillList(x.vec,dims0)
        }else{
           x.vec <- vector("list",dims0)
        }
        if(is(distr, "AbscontDistribution")){
           lower0 <- getLow(distr, eps = getdistrOption("TruncQuantile")*2)
           upper0 <- getUp(distr, eps = getdistrOption("TruncQuantile")*2)
           me <- median(distr); s <- IQR(distr)
           lower1 <- me - 6 * s
           upper1 <- me + 6 * s
           lower <- max(lower0, lower1)
           upper <- min(upper0, upper1)
           h <- upper - lower
           plty <- "l"
           if(is.null(lty)) lty <- "solid"
           for(i in 1:dims0){
               if(is.null(x.vec[[i]])){
                  if(is.finite(xm[i])) lower <- min(lower,xm[i])
                  if(is.finite(xM[i])) upper <- max(upper,xM[i])
                  if(scaleX[i]){
                     xpl <- scaleX.fct[[i]](lower - 0.1*h)
                     xpu <- scaleX.fct[[i]](upper + 0.1*h)
                     xp.vec <- seq(from = xpl, to = xpu, length = 1000)
                     x.vec[[i]] <- scaleX.inv[[i]](xp.vec)
                  }else{
                     x.vec[[i]] <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
                  }
               }
               x.vec[[i]] <- x.vec[[i]][(x.vec[[i]]>=xm[i]) & (x.vec[[i]]<=xM[i])]
           }
        }else{
           for(i in 1:dims0){
              if(!is.null(x.vec[[i]])){
                 if(is(distr, "DiscreteDistribution"))
                 x.vec[[i]] <- intersect(x.vec[[i]],support(distr))
              }else{
                 if(is(distr, "DiscreteDistribution")){
                    x.vec[[i]] <- support(distr)
                 }else{
                    x.vec[[i]] <- r(distr)(1000)
                    x.vec[[i]] <- sort(unique(x.vec[[i]]))
                 }
              }
              x.vec[[i]] <- x.vec[[i]][(x.vec[[i]]>=xm[i]) & (x.vec[[i]]<=xM[i])]
           }
           plty <- "p"
           if(is.null(lty)) lty <- "dotted"
        }
     return(list(x.vec=x.vec, lty = lty, plty = plty))
}

.getXlimYlim <- function(dots,dotsP, dims0, xlim, ylim){
        xm <- rep(-Inf,dims0)
        xM <- rep(Inf,dims0)

        if(!is.null(xlim)){
           if(! length(xlim) %in% c(2,2*dims0))
                stop("Wrong length of Argument xlim");
           xlim <- matrix(xlim, 2,dims0)
           for(i in 1:dims0){
              xm[i] <- min(xlim[,i])
              xM[i] <- max(xlim[,i])
              dotsP[[i]]$xlim <- xlim[,i]
           }
        }

        ym <- rep(-Inf,dims0)
        yM <- rep(Inf,dims0)

        if(!is.null(ylim)){
           if(! length(ylim) %in% c(2,2*dims0))
                stop("Wrong length of Argument ylim");
           ylim <- matrix(ylim, 2,dims0)
           for(i in 1:dims0){
              ym[i] <- min(ylim[,i])
              yM[i] <- max(ylim[,i])
              dotsP[[i]]$ylim <- ylim[,i]
           }
        }
        dots$ylim <- dots$xlim <- NULL
    return(list(dots=dots, dotsP = dotsP, xlim = xlim, ylim = ylim, xm = xm, xM = xM))
}


.prepareTitles <- function(withSubst, presubArg2, presubArg3, dots, mainText,
                           L2Fam, inner, dims0, dims, to.draw, trafO, obj, type, bmar, tmar){

           mainL <- FALSE
           lineT <- NA
           subL <- FALSE
           .mpresubs <- if(withSubst){function(inx)
                    .presubs(inx, presubArg2, presubArg3)
                    } else function(inx)inx

            main <- dots$main
            sub <- dots$sub
            cex.inner <- dots$cex.inner

            if (!is.null(main)){
                 mainL <- TRUE
                 if (is.logical(main)){
                     if (!main) mainL <-  FALSE
                     else
                          main <- mainText ###
                                  ### double  %% as % is special for gettextf
                     }
                 main <- .mpresubs(main)
                 if (mainL) {
                     if(is.null(tmar))
                        tmar <- 5
                     if(is.null(cex.inner))
                        cex.inner <- .65
                     lineT <- 0.6
                     }
             }
             if (hasArg(sub)){
                 subL <- TRUE
                 if (is.logical(sub)){
                     if (!sub) subL <-  FALSE
                     else       sub <- gettextf("generated %%D")
                                  ### double  %% as % is special for gettextf
                 }
                 sub <- .mpresubs(sub)
                 if (subL)
                     if (is.null(bmar)) bmar <- 6
             }
             .mknam <- function(val, iP = "", txt){
                nm <- names(val)
                nms <- if(is.null(nm)) NULL else paste("'", nm, "' = ", sep = "")
                iP <- paste(iP, "\n", gettext(txt), " (",
                        paste(nms, round(val,3), collapse = ", "), ")", sep = "")
             }


             innerParam <- .mknam(L2Fam@param@main, "", "with main parameter")
             if(!is.null(L2Fam@param@nuisance))
                 innerParam <- .mknam(L2Fam@param@nuisance, innerParam,
                                 "and nuisance parameter")
             if(!is.null(L2Fam@param@fixed))
                 innerParam <- .mknam(L2Fam@param@fixed, innerParam,
                             "and fixed known parameter")

             if(!is.logical(inner)){
                #if(!is.character(inner))
                #stop("Argument 'inner' must either be 'logical' or a 'list'")
                if(!is.list(inner))
                    inner <- as.list(inner)
                inner <- lapply(.mpresubs,inner)
                innerT <- .fillList(inner,(type=="info")+dims)
                inf1 <- 0
                if(type=="info"){ if(1 %in% to.draw) inf1 <- 1}
                if(dims0<dims){
                   innerT0 <- innerT
                   for(i in 1:dims0) innerT[inf1+to.draw[i]] <- innerT0[inf1+i]
                }
                innerL <- TRUE
            }else{if(any(is.na(inner))||any(!inner)) {
                     innerT <- as.list(rep("",(type=="info")+dims))
                     innerL <- FALSE
                }else{innerL <- TRUE
                      tnm  <- rownames(trafO)
                      tnms <- if(is.null(tnm)) paste(1:dims) else
                                               paste("'", tnm, "'", sep = "")
                      if(type == "info")
                        innerT <- c( paste(gettext("Absolute information of (partial) IC for "),
                                       name(L2Fam)[1], sep =""),
                                   paste(gettext("Relative information of \ncomponent "),
                                       tnms,
                                       gettext(" of (partial) IC\nfor "),
                                       name(L2Fam)[1], sep =""))
                      if(type=="compare")
                        innerT <- paste(gettext("Component "),  tnms,
                                   gettext(" of (partial) IC\nfor "),
                                   name(L2Fam)[1], sep ="")
                      if(type=="all")
                        innerT <- paste(gettextf("Component "),  tnms,
                                  gettextf("\nof"), #gettextf(" of L_2 derivative\nof"),
                                  name(obj)[1])
                      innerT <- as.list(paste(innerT, innerParam))
                }

           }
  return(list(dots = dots, main = main, mainL = mainL, lineT = lineT,
                           sub = sub, subL = subL, bmar = bmar, tmar = tmar,
                           innerT = innerT, innerL = innerL, .mpresubs = .mpresubs
                           ))
}

.getToDraw <- function(dims, trafO, L2Fam, to.draw.arg, Abs=NULL){
        to.draw <- 1:(dims+!is.null(Abs))
        dimnms  <- c(rownames(trafO))
        if(is.null(dimnms))
           dimnms <- names(main(L2Fam@param))# paste("dim",1:dims,sep="")
        if(is.null(dimnms))
           dimnms <- paste("dim",1:dims,sep="")
        dimnms <- c(Abs,dimnms)
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg))
                 to.draw <- pmatch(to.draw.arg, dimnms)
            else if(is.numeric(to.draw.arg))
                 to.draw <- to.draw.arg
        }
   return(to.draw)
}

.preparePanelFirstLast <- function(with.automatic.grid , dims0, pF.0, pL.0,
    logArg, scaleX, scaleY, x.ticks, y.ticks, scaleX.fct, scaleY.fct){

            with.automatic.grid <- rep(with.automatic.grid, length.out=dims0)

            pF <- pF.L <- .fillList(pF.0,dims0)
            if(dims0) for(i in 1:dims0){
                     if(!is.null(logArg)){
                        scaleX[i] <- scaleX[i] & !grepl("x",logArg[i])
                        scaleY[i] <- scaleY[i] & !grepl("y",logArg[i])
                     }
                     if(with.automatic.grid[i]&&
                        (scaleX[i]||scaleY[i])
                     ){
                        pF[[i]] <-  .producePanelFirstSn(
                             pF.L[[i]],
                              x.ticks = x.ticks[[i]], scaleX = scaleX[i],
                                  scaleX.fct = scaleX.fct[[i]],
                              y.ticks = y.ticks[[i]], scaleY = scaleY[i],
                                  scaleY.fct = scaleY.fct[[i]])
                     }
            }

            gridS <- if(any(with.automatic.grid))
                 substitute({grid <- function(...){}}) else expression({})

            pL <- pL.0
            if(is.list(pL.0)){
               pL.0 <- .fillList(pL.0, dims0)
               if(dims0) for(i in 1:dims0){
                    if(is.null(pL.0[[i]])) pL.0[[i]] <- expression({})
               }
               pL <- substitute({pL1 <- pL0
                           pL1[[i]]},
                           list(pL0=pL.0))
            }
            return(list(pF=pF, pL=pL, gridS = gridS))
}
