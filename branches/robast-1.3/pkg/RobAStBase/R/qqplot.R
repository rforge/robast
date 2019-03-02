################################################################
# QQ - Plot functions in package RobAStBase
################################################################



.fadeColor <- function(col,x, bg = "white"){
 ind <- seq(along=x)
 col <- .makeLenAndOrder(col,ind)
 colx <- t(sapply(ind,function(i) colorRamp(c(bg,col[i]))(x[i])))
 colv2col <- function(colvec)
   rgb(red = colvec[1], green = colvec[2], blue = colvec[3], maxColorValue = 255)
 apply(colx,1,function(x) colv2col(x))
}

## into RobAStBase
setMethod("qqplot", signature(x = "ANY",
                              y = "RobModel"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ..., distance = NormType(),
    n.adj = TRUE){

    mc <- match.call(call = sys.call(sys.parent(1)))
    dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
    args0 <- list(x = x, y = y,
                  n = if(!missing(n)) n else length(x),
                  withIdLine = withIdLine, withConf = withConf,
    withConf.pw  = if(!missing(withConf.pw)) withConf.pw else if(!missing(withConf)) withConf else NULL,
    withConf.sim = if(!missing(withConf.sim)) withConf.sim else if(!missing(withConf)) withConf else NULL,
                  plot.it = plot.it, xlab = xlab, ylab = ylab, distance=distance, n.adj=n.adj)

    plotInfo <- list(call=mc, dots=dots, args=args0)

    xcc <- as.character(deparse(mc$x))
    if(missing(xlab)) mc$xlab <- xcc
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]

    if(is.null(mcl$n.CI)) mcl$n.CI <- n
    if(n.adj){
       r <- radius(neighbor(y))
       n <- floor((1-r)*n)
    }
    if(is.null(mcl$alpha.CI))
       mcl$alpha.CI <- .95
    cor <- radius(neighbor(y))
    mcl$legend.alpha <- eval(mcl$alpha.CI)
    mcl$alpha.CI <- min(eval(mcl$alpha.CI)+cor,1)


    mcl$n <- n
    mcl$y <- y@center
    mcl$legend.pref <- paste(mcl$legend.pref,"outlier-adjusted",sep="")


    xD <- fct(distance)(x)
    x.cex <- 3/(1+log(1+xD))
    mcl$cex.pts <- x.cex

    retv <- do.call(getMethod("qqplot", signature(x="ANY", y="ProbFamily")),
            args=mcl)
    retv$call <- retv$dots <- retv$args <- NULL
    plotInfo <- c(plotInfo,retv)
    class(plotInfo) <- c("qqplotInfo","DiagnInfo")
    return(invisible(plotInfo))
    })


## into RobAStBase
setMethod("qqplot", signature(x = "ANY", y = "InfRobModel"),
      function(x, y, n = length(x), withIdLine = TRUE, withConf = TRUE,
               withConf.pw  = withConf,  withConf.sim = withConf,
               plot.it = TRUE, xlab = deparse(substitute(x)),
               ylab = deparse(substitute(y)), ..., cex.pts.fun = NULL, n.adj = TRUE){

    mc <- match.call(call = sys.call(sys.parent(1)))
    dots <- match.call(call = sys.call(sys.parent(1)),
                       expand.dots = FALSE)$"..."
    args0 <- list(x = x, y = y,
                  n = if(!missing(n)) n else length(x),
                  withIdLine = withIdLine, withConf = withConf,
    withConf.pw  = if(!missing(withConf.pw)) withConf.pw else if(!missing(withConf)) withConf else NULL,
    withConf.sim = if(!missing(withConf.sim)) withConf.sim else if(!missing(withConf)) withConf else NULL,
                  plot.it = plot.it, xlab = xlab, ylab = ylab, cex.pts.fun=cex.pts.fun,
                  n.adj = n.adj)

    plotInfo <- list(call=mc, dots=dots, args=args0)
    if(missing(xlab)) mc$xlab <- as.character(deparse(mc$x))
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]
    if(is.null(mcl$distance)) distance <- NormType()

    if(is.null(mcl$alpha.CI))
       mcl$alpha.CI <- .95
    cor <- radius(neighbor(y))/sqrt(n)
    mcl$legend.alpha <- eval(mcl$alpha.CI)
    mcl$alpha.CI <- min(eval(mcl$alpha.CI)+cor,1)



    mcl$n <- n

    if(is.null(mcl$n.CI)) mcl$n.CI <- n
    if(n.adj){
       r <- radius(neighbor(y))
       mcl$n.CI <- floor((1-r/sqrt(n))*n)
    }
    mcl$y <- y@center
    mcl$legend.pref <- paste(mcl$legend.pref,"outlier-adjusted",sep="")
    
    FI <- PosSemDefSymmMatrix(FisherInfo(y@center))
    L2D <- as(diag(nrow(FI)) %*% L2deriv(y@center), "EuclRandVariable")
    L2Dx <- evalRandVar(L2D,matrix(x))[,,1]
    scx <-  distr::solve(sqrt(FI),L2Dx)
    xD <- fct(distance)(scx)
    cex.pts <- if(is.null(mcl[["cex.pts"]])){
                  if(is.null(mcl[["cex"]])){
                     par("cex")
                  }else{
                     eval(mcl$cex)}
               }else{
                  eval(mcl$cex.pts)
               }

    x.cex <- 3/(1+.cexscale(xD,xD,cex=cex.pts, fun = cex.pts.fun))

    mcl$cex.pts <- x.cex

    retv <- do.call(getMethod("qqplot", signature(x="ANY", y="ProbFamily")),
            args=mcl)
    retv$call <- retv$dots <- retv$args <- NULL
    plotInfo <- c(plotInfo,retv)
    class(plotInfo) <- c("qqplotInfo","DiagnInfo")
    return(invisible(plotInfo))
    })

## into RobAStBase
setMethod("qqplot", signature(x = "ANY",
                              y = "kStepEstimate"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...,
    exp.cex2.lbs = -.15,
    exp.cex2.pts = -.35,
    exp.fadcol.lbs = 1.85,
    exp.fadcol.pts = 1.85,
    bg = "white"
    ){

    args0 <- list(x=x,y=y,n=n,withIdLine=withIdLine, withConf=withConf,
        withConf.pw  = if(!missing(withConf.pw)) withConf.pw else if(!missing(withConf)) withConf else NULL,
        withConf.sim = if(!missing(withConf.sim)) withConf.sim else if(!missing(withConf)) withConf else NULL,
        plot.it = plot.it, xlab = xlab, ylab = ylab, exp.cex2.lbs=exp.cex2.lbs,
        exp.cex2.pts=exp.cex2.pts, exp.fadcol.lbs=exp.fadcol.lbs,
        exp.fadcol.pts=exp.fadcol.pts, bg=bg)

    mc <- match.call(call = sys.call(sys.parent(1)))
    mc1 <- match.call(call = sys.call(sys.parent(1)), expand.dots=FALSE)
    dots <- mc1$"..."
    plotInfo <- list(call=mc, dots=dots, args=args0)

    if(missing(xlab)) mc$xlab <- as.character(deparse(mc$x))
    if(missing(ylab)) mc$ylab <- as.character(deparse(mc$y))
    mcl <- as.list(mc)[-1]

    IC <- pIC(y)
    if(!is(IC,"IC"))
       stop("IC of the kStepEstimator needs to be of class 'IC'")

    L2Fam <- eval(IC@CallL2Fam)
    param <- ParamFamParameter(main=untransformed.estimate(y), nuisance=nuisance(y),
                               fixed=fixed(y))
    L2Fam0 <- modifyModel(L2Fam,param)
    mcl$y <- L2Fam0

    if(is(IC,"HampIC")){
      dim0 <- nrow(FisherInfo(L2Fam))
      L <- as(diag(dim0)%*%L2Fam@L2deriv, "EuclRandVariable")
      L.fct <- function(x) evalRandVar(L,x)

      w.fct <- function(x)
               weight(weight(IC))(L.fct(matrix(x))[,,1])

      wx <- 1/(1+w.fct(x))
      if(max(wx)>1) wx <- wx/max(wx)
      mcl$order.traf <- function(x) 1/w.fct(x)

      cex.lbs <- if(is.null(mcl$cex.lbs))  par("cex")  else eval(mcl$cex.lbs)
      cex.pts <- if(is.null(mcl$cex.pts))  par("cex")  else eval(mcl$cex.pts)
      mcl$cex.lbs <- cex.lbs*wx^exp.cex2.lbs
      mcl$cex.pts <- cex.pts*wx^exp.cex2.pts

      col.lbs <- if(is.null(mcl$col.lbs))  par("col")  else eval(mcl$col.lbs)
      col.pts <- if(is.null(mcl$col.pts))  par("col")  else eval(mcl$col.pts)
      mcl$col.lbs <- .fadeColor(col.lbs,wx^exp.fadcol.lbs, bg = bg)
      mcl$col.pts <- .fadeColor(col.pts,wx^exp.fadcol.pts, bg = bg)
    }

    retv <- do.call(getMethod("qqplot", signature(x="ANY", y="ProbFamily")),
            args=mcl)
    retv$call <- retv$dots <- retv$args <- NULL
    plotInfo <- c(plotInfo,retv)
    class(plotInfo) <- c("qqplotInfo","DiagnInfo")
    return(invisible(plotInfo))
    })
