setMethod("comparePlot", signature("IC","IC"),
    function(obj1,obj2, obj3 = NULL, obj4 = NULL, data = NULL,
             ..., withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col = par("col"), lwd = par("lwd"), lty, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], 
             with.legend = TRUE, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, col.pts = par("col"),
             pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
             lab.pts = NULL, lab.font = NULL,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE){

        mcl <- match.call(call = sys.call(sys.parent(1)))

        L2Fam <- eval(x@CallL2Fam); trafO <- trafo(L2Fam@param)
        dims  <- nrow(trafO); to.draw <- 1:dims
        if(! is.null(to.draw.arg)){
            if(is.character(to.draw.arg))
                 to.draw <- pmatch(to.draw.arg, dimnms)
            else if(is.numeric(to.draw.arg))
                 to.draw <- to.draw.arg
        }
        dims0 <- length(to.draw)

        MBRB <- matrix(rep(t(MBRB), length.out=dims0*2),ncol=2, byrow=T)
        if(withMBR && all(is.na(MBRB))){
           robModel <- InfRobModel(center = L2fam, neighbor =
                             ContNeighborhood(radius = 0.5))
           ICmbr <- try(optIC(model = robModel, risk = asBias()), silent=TRUE)
           if(!is(ICmbr,"try-error"))
              MBRB <- .getExtremeCoordIC(ICmbr, distribution(L2Fam), todraw)
           else withMBR <- FALSE
        }
        mcl$MBRB <- MBRB
        mcl$withMBR <- withMBR
        do.call(getMethod("comparePlot", signature("IC","IC"),
                           where="RobAStBase"), mcl)
    })
