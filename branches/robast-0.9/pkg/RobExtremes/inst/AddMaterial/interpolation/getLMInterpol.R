getLMs <- function(Gridnam,Famnam,xi=0.7, baseDir="C:/rtest/robast"){
   ## Gridnam in (Sn,OMSE,RMXE,MBRE)
   ## Famnam in "Generalized Pareto Family",
   ##           "Generalized Extreme Value Family with positive shape parameter: Frechet Family",
   ##           "Gamma family",
   ##           "Weibull Family"
   ## xi Scaleparameter (can be vector)
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
   fct <- get(GN,envir=nE)[[Famnam]]$fct

   if(!isSn)){
   ## für Gridnam != Sn ist LM für jeden xi Wert ein Vektor der Länge 13, genauer
   #           in 1:13 (clip=b, cent.a=a1.a,a2.a, cent.i=a1.i,a2.i,
   ##                  stand.a=A.a=matrix(c(A11.a,(A12.a+A21.a)/2,
   #                                       (A12.a+A21.a)/2,A.22.a), 2, 2),
   ##                  stand.i=A.i=matrix(c(A11.i,(A12.i+A21.i)/2,
   #                                       (A12.i+A21.i)/2,A.22.i), 2, 2),
   ##                 und optIC = Y.a min(1,b/norm(Y.i)), Y.* = A.* Lambda - a.*
      len <- length(fct)
      LM <- sapply(1:len, function(i) fct[[i]](xi))
      if(length(xi)==1) LM <- matrix(LM,ncol=len)
      colnames(LM) <- c("b","a1.a", "a2.a", "a1.i", "a2.i", "A11.a",
                 "A12.a", "A21.a", "A22.a", "A11.i", "A12.i", "A21.i", "A22.i")
      return(cbind(xi,LM))
   }else{
      Sn <- fct(xi)
      return(cbind(xi,Sn))
   }
}
