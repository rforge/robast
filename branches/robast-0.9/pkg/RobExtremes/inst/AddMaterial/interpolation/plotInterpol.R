plotLM <- function(Gridnam,Famnam,whichLM, baseDir="C:/rtest/robast",
               withSmooth = FALSE, plotGridRestriction = NULL,
               smoothtry = FALSE, df = NULL, gridRestrForSmooth = NULL,
               prehook={}, posthook={}, ...){
   ## Gridnam in (Sn,OMSE,RMXE,MBRE) ## uses partial matching!!
   ## Famnam in "Generalized Pareto Family", ## uses partial matching!!
   ##           "GEV Family",
   ##           "Gamma family",
   ##           "Weibull Family"
   ## whichLM  is ignored for Gridnam == Sn
   #           in 1:13 (clip=b, cent.a=a1.a,a2.a, cent.i=a1.i,a2.i,
   ##                  stand.a=A.a=matrix(c(A11.a,(A12.a+A21.a)/2,
   #                                       (A12.a+A21.a)/2,A.22.a), 2, 2),
   ##                  stand.i=A.i=matrix(c(A11.i,(A12.i+A21.i)/2,
   #                                       (A12.i+A21.i)/2,A.22.i), 2, 2),
   ##                 and optIC = Y.a min(1,b/norm(Y.i)), Y.* = A.* Lambda - a.*
   ##          or "all" then all LMs are plotted
   ## basedir: folder with r-forge svn checkout
   ## plotGridRestriction: an expression that can be used as index in
   ##                      xi[plotGridRestriction] to restrict the plotted
   ##                      grid-values
   ## prehook: an expression to be evaluated before plotting --- typically something
   ##          like pdf("myfile.pdf")
   ## posthook: an expression to be evaluated after plotting --- typically something
   ##          like dev.off()
   ## withSmooth: logical; shall item grid or gridS be used for plotting
   ## ---------------------------
   ###  for interactive try-out of several smoothing values
   ##
   ## smoothtry: logical; shall interactive try-out of smoothing be used
   ##            if TRUE overrides withSmooth
   ## df: smoothing parameter (see below)
   ## gridRestrForSmooth: restriction of smoothing for particular theta-grid-values
   ##        (see below)
   ###
   ## copied from help to .MakeSmoothGridList
   ##
#   \item{df}{argument \code{df} of \code{\link{smooth.spline}}; if \code{NULL}
#            (default) it is omitted (and the default of
#            \code{\link{smooth.spline}} used); controls the degree to which
#            we smooth; can be vectorized; to allow for \code{NULL}-entries
#            in some (of the 13) LMs, it can also be a list of length 13,
#            some entries being \code{NULL}, some numeric. }
#  \item{gridRestrForSmooth}{an expression that can be used as index in
#     \code{theta[gridRestrForSmooth]} to restrict the grid-values to
#     be smoothed; the excluded grid values are left unchanged. If the argument
#     is \code{NULL} no restriction is used. Can be a matrix of same dimension
#     as the \code{Y}-grid to allow for column-individual restrictions,
#     or a list of same length as number of columns of \code{Y}
#     with columnwise restrictions of \code{Y} (and \code{NULL} entries
#     are interpreted as no restriction). }

   file <- file.path(baseDir, "branches/robast-0.9/pkg/RobAStRDA/R/sysdata.rda")
   if(!file.exists(file)) stop("Fehler mit Checkout")
   nE <- new.env()
   load(file, envir=nE)
   Gnams <- c("Sn","OMSE","RMXE","MBRE")
   Fnams <- c("Generalized Pareto Family",
              "GEV Family",
              "Gamma family",
              "Weibull Family")
   Gridnam <- Gnams[pmatch(Gridnam, Gnams)]
   Famnam <- Fnams[pmatch(Famnam, Fnams)]
   if(! Gridnam %in% Gnams) stop("Falscher Gittername")
   if(! Famnam %in% Fnams) stop("Falscher Familienname")
   isSn <- (Gridnam == "Sn")
   Famnam0 <- gsub(" ","",Famnam)
   GN0 <- Gridnam; if(isSn) GN0 <- "SnGrids"
   GN <- paste(".",GN0,sep="")
   funN <- paste("fun",".",if(getRversion()<"2.16") "O" else "N",sep="")
   if(!smoothtry){
      gN <- if(withSmooth) "gridS" else "grid"
      gr <- get(GN,envir=nE)[[Famnam0]][[gN]]
   }else{
     gr <- get(GN,envir=nE)[[Famnam0]][["grid"]]
#     gr <- RobAStRDA:::.MakeSmoothGridList(gr[,1],gr[,-1], df = df,
#                        gridRestrForSmooth = gridRestrForSmooth)
     gr <- .MakeSmoothGridList(gr[,1],gr[,-1], df = df,
                        gridRestrForSmooth = gridRestrForSmooth)
   }
   if(is.null(plotGridRestriction)) plotGridRestriction <- rep(TRUE, nrow(gr))
   if(!isSn) if(whichLM!="all") if(whichLM<1 | whichLM>13) stop("Falsche Koordinate")
   if(!isSn) if(whichLM=="all"){
      eval(prehook)
      par(mfrow=c(4,4))
      for(i in 2:14)
          plot(gr[plotGridRestriction,1], gr[plotGridRestriction,i], ...)
      par(mfrow=c(1,1))
      eval(posthook)
   return(invisible(NULL))
   }
   if(isSn) whichLM <- 1
   wM <- whichLM + 1
   eval(prehook)
   plot(gr[gridRestriction,1], gr[gridRestriction,wM], ...)
   eval(posthook)
   return(invisible(NULL))
}

if(FALSE){
## Examples
plotLM("OMSE","Gamma","all", type="l", plotG=-(1:8), main ="Gamma-OMSE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("OMSE","Gener","all", type="l", plotG=-(1:8), main ="GPD-OMSE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("OMSE","GEV","all", type="l", plotG=-(1:8), main ="GEV-OMSE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("OMSE","Wei","all", type="l", plotG=-(1:8), main ="Weibull-OMSE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("MBRE","Gam","all", type="l", plotG=-(1:8), main ="Gamma-MBRE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("MBRE","Gene","all", type="l", plotG=-(1:8), main ="GPD-MBRE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("MBRE","GE","all", type="l", plotG=-(1:8), main ="GEV-MBRE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("MBRE","Wei","all", type="l", plotG=-(1:8), main ="Weibull-MBRE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("RMXE","Gam","all", type="l", plotG=-(1:8), main ="Gamma-RMXE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("RMXE","Gene","all", type="l", plotG=-(1:8), main ="GPD-RMXE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("RMXE","GE","all", type="l", plotG=-(1:8), main ="GEV-RMXE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("RMXE","Wei","all", type="l", plotG=-(1:8), main ="Weibull-RMXE", xlab=expression(xi),ylab="LM", pre=windows())
plotLM("MBRE","GE","all", type="l")
plotLM("MBRE","GE","all", type="l", sm = TRUE, df = 10, gridR = -(1:15))
plotLM("MBRE","GE","all", type="l", sm = TRUE, df = 4, gridR = -(1:15))
}