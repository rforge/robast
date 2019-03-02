.MakeSmoothGridList <- RobExtremes:::.MakeSmoothGridList
plotLM <- function(Gridnam,Famnam,whichLM, baseDir="C:/rtest/robast",
               withSmooth = FALSE, plotGridRestriction = NULL,
               smoothtry = FALSE, df = NULL, gridRestrForSmooth = NULL,
               prehook={}, posthook={}, ylab=NULL, xlab=NULL, main = NULL,
               lwd=NULL, lty= NULL, col =NULL, inputSmooth = FALSE, ...,
               rdaFilen = "sysdata.rda",
               rdaRelPath = "branches/robast-1.1/pkg/RobExtremesBuffer"){
   ## Gridnam in (Sn,OMSE,RMXE,MBRE) ## uses partial matching!!
   ## Famnam in "Generalized Pareto Family", ## uses partial matching!!
   ##           "GEV Family",
   ##           "GEVU Family",
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
   ## withSmooth: logical; shall item grid or gridS be used for plotting / smoothing
   ##                      this way smoothing can be split up in several steps
   ## ---------------------------
   ###  for interactive try-out of several smoothing values
   ##
   ## smoothtry: logical; shall interactive try-out of smoothing be used
   ##            if TRUE overrides withSmooth
   ## df: smoothing parameter (see below)
   ## gridRestrForSmooth: restriction of smoothing for particular theta-grid-values
   ##        (see below)
   ## ylab, xlab, lty, lwd, col parameters for plot (or NULL, then defaults are used)
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

   file <- file.path(baseDir, rdaRelPath, rdaFilen)
#   file <- file.path(baseDir, "branches/robast-1.0/pkg/RobAStRDA/R/sysdata.rda")
   if(!file.exists(file)) stop("Fehler mit Checkout")
   nE <- new.env()
   load(file, envir=nE)
   Gnams <- c("Sn","OMSE","RMXE","MBRE")
   Fnams <- c("Generalized Pareto Family",
              "GEV Family",
              "GEVU Family",
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
   ## gN <- if(withSmooth) "gridS" else "grid"
   gr0 <- get(GN,envir=nE)[[Famnam0]][["grid"]]
   gr1 <- get(GN,envir=nE)[[Famnam0]][["gridS"]]
   gr <- if(withSmooth) gr1 else gr0
   if(smoothtry){
     gr <- .MakeSmoothGridList(gr[,1],gr[,-1], df = df,
                        gridRestrForSmooth = gridRestrForSmooth)
   }else{gr <- gr1}
   print(round(head(gr0[,1]),3))
   print(round(tail(gr0[,1]),3))
   print(round((gr0[gr0[,1]>5,1]),3))
   print(round(summary(diff(gr0[,1])),3))
   print(c("n"=sum(!is.na(gr0[,1])),"NA"=sum(is.na(gr0[,1]))))
   z0 <- if(Famnam=="GEVU Family") 25 else 13
   z <- if(whichLM=="all")  z0 else 1
   if(is.null(plotGridRestriction)){
      plotGridRestriction <- list(rep(TRUE, nrow(gr)))
      pl<- vector("list",z)
      if(z>1) pl[1:z] <- plotGridRestriction else pl[1] <- plotGridRestriction
      plotGridRestriction <- pl
   }else{
      pl<- vector("list", z)
      pla <- if(is.list(plotGridRestriction)) plotGridRestriction else list(plotGridRestriction)
      pl1 <- rep(pla, length.out=z)
      plotGridRestriction <- pl1
   }
   namesLM <- c("b","a.a[sc]","a.a[sh]","z.i[sc]","z.i[sh]",
                "A.a[sc,sc]","A.a[sc,sh]","A.a[sh,sc]","A.a[sh,sh]",
                "A.i[sc,sc]","A.i[sc,sh]","A.i[sh,sc]","A.i[sh,sh]")
   if(!isSn) if(whichLM!="all") if(whichLM<1 | whichLM>z0) stop("Falsche Koordinate")
   if(missing(ylab)||is.null(ylab)) ylab <- "LM"
   if(missing(xlab)||is.null(xlab)) xlab <- expression(xi)
   if(missing(main)||is.null(main)) main <- paste(Gridnam,gsub(" [F,f]amily","",Famnam),sep="-")
   if(missing(lty)||is.null(lty)) lty <- c(2,3,1)
   if(missing(lwd)||is.null(lwd)) lwd <- c(0.8,0.8,1.8)
   if(missing(col)||is.null(col)) col <- 1:3
   if(!isSn) if(whichLM=="all"){
      eval(prehook)
      if(z0==25) par(mfrow=c(6,6)) else par(mfrow=c(4,4))    
      if(z0==25) upP=26 else upP=14
      for(i in 2:upP){
          pla <- plotGridRestriction[[i-1]]
          if(is.null(pla)) pla <- 1:nrow(gr)
          matplot(gr[pla,1], cbind(gr0[pla,i],gr[pla,i]), xlab=xlab, type="n",
                  ylab=paste(ylab,namesLM[i-1]), main=main,  ...)
          matlines(gr[pla,1],
             cbind(gr0[pla,i],gr1[pla,i],gr[pla,i]),lwd=lwd, lty=lty, col=col)
      }
      par(mfrow=c(1,1))
      eval(posthook)
   return(invisible(NULL))
   }
   if(isSn) whichLM <- 1
   wM <- whichLM + 1
   eval(prehook)
   pla <- plotGridRestriction[[1]]
   if(is.null(pla)) pla <- 1: nrow(gr)
   print(wM)
   print(head(gr[gridRestrForSmooth[[1]],1]))
   matplot(gr[pla,1], cbind(gr0[pla,wM],gr[pla,wM]), type="n",
            xlab=xlab, ylab=paste(ylab,namesLM[wM-1]), main = main, ...)
   matlines(gr[pla,1],
             cbind(gr0[pla,wM],gr1[pla,wM],gr[pla,wM]),lwd=lwd, lty=lty, col=col)
   eval(posthook)
   return(invisible(NULL))
}

if(FALSE){
## Examples
plotLM("OMSE","Gamma","all", plotG=-(1:8), pre=windows())
plotLM("OMSE","Gener","all", plotG=-(1:8), pre=windows())
plotLM("OMSE","GEV","all", plotG=-(1:8), pre=windows())
plotLM("OMSE","Wei","all", plotG=-(1:8), pre=windows())
plotLM("MBRE","Gam","all", plotG=-(1:8), pre=windows())
plotLM("MBRE","Gene","all", plotG=-(1:8), pre=windows())
plotLM("MBRE","GE","all", plotG=-(1:8), pre=windows())
plotLM("MBRE","Wei","all", plotG=-(1:8), pre=windows())
plotLM("RMXE","Gam","all", plotG=-(1:8), pre=windows())
plotLM("RMXE","Gene","all", plotG=-(1:8), pre=windows())
plotLM("RMXE","GE","all", plotG=-(1:8), pre=windows())
plotLM("RMXE","Wei","all", plotG=-(1:8), pre=windows())
plotLM("MBRE","GE","all", type="l")
plotLM("MBRE","GE","all", sm = TRUE, df = 10, gridR = -(1:15))
plotLM("MBRE","GE","all", sm = TRUE, df = 4, gridR = -(1:15))
plotLM("OMSE","Gamma","all", plotG=-(1:8), withS=TRUE, pre=pdf("Gamma-OMSE-s.pdf"), post=dev.off())
plotLM("OMSE","Gener","all", plotG=-(1:8), withS=TRUE, pre=pdf("GPD-OMSE-s.pdf"), post=dev.off())
plotLM("OMSE","GEV","all", plotG=-(1:8), withS=TRUE, pre=pdf("GEV-OMSE-s.pdf"), post=dev.off())
plotLM("OMSE","Wei","all", plotG=-(1:8), withS=TRUE, pre=pdf("Weibull-OMSE-s.pdf"), post=dev.off())
plotLM("MBRE","Gam","all", plotG=-(1:8), withS=TRUE, pre=pdf("Gamma-MBRE-s.pdf"), post=dev.off())
plotLM("MBRE","Gene","all", plotG=-(1:8), withS=TRUE, pre=pdf("GPD-MBRE-s.pdf"), post=dev.off())
plotLM("MBRE","GE","all", plotG=-(1:8), withS=TRUE, pre=pdf("GEV-MBRE-s.pdf"), post=dev.off())
plotLM("MBRE","Wei","all", plotG=-(1:8), withS=TRUE, pre=pdf("Weibull-MBRE-s.pdf"), post=dev.off())
plotLM("RMXE","Gam","all", plotG=-(1:8), withS=TRUE, pre=pdf("Gamma-RMXE-s.pdf"), post=dev.off())
plotLM("RMXE","Gene","all", plotG=-(1:8),withS=TRUE,  pre=pdf("GPD-RMXE-s.pdf"), post=dev.off())
plotLM("RMXE","GE","all", plotG=-(1:8), withS=TRUE, pre=pdf("GEV-RMXE-s.pdf"), post=dev.off())
plotLM("RMXE","Wei","all", plotG=-(1:8), withS=TRUE, pre=pdf("Weibull-RMXE-s.pdf"), post=dev.off())
}
