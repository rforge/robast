getLMs <- function(Gridnam,Famnam,xi=0.7, baseDir="C:/rtest/robast", withPrint=FALSE, withLoc = FALSE){
   ## Gridnam in (Sn,OMSE,RMXE,MBRE) ## uses partial matching!!
   ## Famnam in "Generalized Pareto Family",
   ##           "GEV Family",
   ##           "Gamma family",
   ##           "Weibull Family"  ## uses partial matching!!
   ## xi Scaleparameter (can be vector)
   ## basedir: Oberverzeichnis des r-forge svn checkouts
   ## withPrint: diagnostischer Output?
   ## withLoc: anzuschalten bei GEVFamilyMuUnknown...
   file <- file.path(baseDir, "branches/robast-0.9/pkg/RobAStRDA/R/sysdata.rda")
   if(!file.exists(file)) stop("Fehler mit Checkout")
   nE <- new.env()
   load(file, envir=nE)
   Gnams <- c("Sn","OMSE","RMXE","MBRE")
   Fnams <- c("Generalized Pareto Family",
              "GEVU Family",
              "GEV Family",
              "Gamma family",
              "Weibull Family")
   Gridnam <- Gnams[pmatch(Gridnam, Gnams)]
   Famnam <- Fnams[pmatch(Famnam, Fnams)]
   if(! Gridnam %in% Gnams) stop("Falscher Gittername")
   if(! Famnam %in% Fnams) stop("Falscher Familienname")
   Famnam0 <- gsub(" ","",Famnam)
   isSn <- (Gridnam == "Sn")
   GN0 <- Gridnam; if(isSn) GN0 <- "SnGrids"
   GN <- paste(".",GN0, sep="")
   funN <- paste("fun",".",if(getRversion()<"2.16") "O" else "N",sep="")
   if(withPrint) print(c(GN, Famnam0, funN))
   fct <- get(GN,envir=nE)[[Famnam0]][[funN]]

   if(!isSn){
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
      if(withLoc){
         colnames(LM) <- c("b","a1.a", "a2.a", "a3.a", "a1.i", "a2.i", "a3.i",
                           "A11.a", "A12.a", "A13.a", "A21.a", "A22.a", "A23.a",
                           "A31.a", "A32.a", "A33.a", "A11.i", "A12.i", "A13.i",
                           "A21.i", "A22.i", "A23.i", "A31.i", "A32.i", "A33.i")
      }else{
         colnames(LM) <- c("b","a1.a", "a2.a", "a1.i", "a2.i", "A11.a",
                           "A12.a", "A21.a", "A22.a", "A11.i", "A12.i",
                           "A21.i", "A22.i")
      }
      return(cbind(xi,LM))
   }else{
      Sn <- fct(xi)
      return(cbind(xi,Sn))
   }
}
