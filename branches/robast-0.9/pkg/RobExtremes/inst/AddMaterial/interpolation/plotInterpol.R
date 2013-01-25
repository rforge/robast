plotLM <- function(Gridnam,Famnam,whichLM, baseDir="C:/rtest/robast"){
   ## Gridnam in (Sn,OMSE,RMXE,MBRE)
   ## Famnam in "Generalized Pareto Family",
   ##           "Generalized Extreme Value Family with positive shape parameter: Frechet Family",
   ##           "Gamma family",
   ##           "Weibull Family"
   ## whichLM  ignoriert für Gridnam == Sn
   #           in 1:13 (clip=b, cent.a=a1.a,a2.a, cent.i=a1.i,a2.i,
   ##                  stand.a=A.a=matrix(c(A11.a,(A12.a+A21.a)/2,
   #                                       (A12.a+A21.a)/2,A.22.a), 2, 2),
   ##                  stand.i=A.i=matrix(c(A11.i,(A12.i+A21.i)/2,
   #                                       (A12.i+A21.i)/2,A.22.i), 2, 2),
   ##                 und optIC = Y.a min(1,b/norm(Y.i)), Y.* = A.* Lambda - a.*
   ## basedir: Oberverzeichnis des r-forge svn checkouts
   file <- file.path(baseDir, "branches/robast-0.9/pkg/RobExtremes/R/sysdata.rda")
   if(!file.exists(file)) stop("Fehler mit Checkout")
   nE <- new.env()
   load(file, envir=nE)
   Gnams <- c("Sn","OMSE","RMXE","MBRE")
   Fnams <- c("Generalized Pareto Family",
              "Generalized Extreme Value Family with positive shape parameter: Frechet Family",
              "Gamma family",
              "Weibull Family")
   if(! Gridnam %in% Gnams) stop("Falscher Gittername")
   if(! Famnam %in% Fnams) stop("Falscher Familienname")
   isSn <- (Gridnam == "Sn")
   GN0 <- Gridnam; if(isSn) GN0 <- "SnGrids"
   GN <- paste(".",GN0,".",if(getRversion()<"2.16") "O" else "N", sep="")
   gr <- get(GN,envir=nE)[[Famnam]]$grid

   if(!isSn) if(whichLM<1 | whichLM>13) stop("Falsche Koordinate")
   if(isSn) whichLM <- 1
   wM <- whichLM + 1
   plot(gr[,1], gr[,wm])
}
