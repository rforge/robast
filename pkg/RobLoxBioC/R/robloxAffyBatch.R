###############################################################################
## Use robloxbioc to preprocess Affymetrix-data - comparable to MAS 5.0
###############################################################################
setMethod("robloxbioc", signature(x = "AffyBatch"),
    function(x, pmcorrect = "roblox", verbose = TRUE,
            eps = NULL, eps.lower = 0, eps.upper = 0.1, steps = 1L, mad0 = 1e-4,
            contrast.tau = 0.03, scale.tau = 10, delta = 2^(-20)) {
        n <- length(x)
        ids <- featureNames(x)
        m <- length(ids)

        CDFINFO <- getCdfInfo(x)
        INDEX <- sapply(ids, get, envir = CDFINFO)
        NROW <- unlist(lapply(INDEX, nrow))
        nr <- as.integer(names(table(NROW)))

        diff.log2 <- function(INDEX, x){
            l.pm <- INDEX[,1]
            if(ncol(INDEX) == 2)
                l.mm <- INDEX[,2]
            else
                l.mm <- integer()
            log2(x[l.pm, , drop = FALSE]) - log2(x[l.mm, , drop = FALSE])
        }
        pm.only <- function(INDEX, x){
            l.pm <- INDEX[,1]
            x[l.pm, , drop = FALSE]
        }

        intensData <- intensity(x)
        rob.est <- matrix(NA, ncol = 2, nrow = m*n)
        if(verbose) cat("PM/MM correcting ...")
        if(pmcorrect == "roblox"){
            res <- lapply(INDEX, diff.log2, x = intensData)
            rob.est1 <- matrix(NA, ncol = 2, nrow = m*n)
            for(k in nr){
                ind <- which(NROW == k)
                temp <- matrix(do.call(rbind, res[ind]), nrow = k)
                ind1 <-  as.vector(sapply(seq_len(n)-1, function(x, ind, m){ ind + x*m }, ind = ind, m = m))
                rob.est1[ind1, 1:2] <- robloxbioc(t(temp), eps = eps, eps.lower = eps.lower, eps.upper = eps.upper, 
                                                  steps = steps, mad0 = mad0)
            }
            sb <- matrix(rob.est1[,1], nrow = m)
            for(k in seq_len(m)){
                IDX <- INDEX[[k]]
                l.pm <- IDX[,1]
                if(ncol(IDX) == 2){
                    l.mm <- IDX[,2]
                }else{
                    l.mm <- integer()
                }
                pps.pm <- intensData[l.pm, , drop = FALSE]
                pps.mm <- intensData[l.mm, , drop = FALSE]
                pps.im <- pps.mm
                l <- t(t(pps.mm >= pps.pm) & (sb[k,] > contrast.tau))
                pps.im[l] <- t(t(pps.pm)/2^sb[k,])[l]
                l <- t(t(pps.mm >= pps.pm) & (sb[k,] <= contrast.tau))
                pps.im[l] <- t(t(pps.pm)/2^(contrast.tau/(1 + (contrast.tau - sb[k,])/scale.tau)))[l]
                pm.corrected <- pmax.int(pps.pm - pps.im, delta)
    #            pm.corrected[pm.corrected < delta] <- delta
                res[[k]] <- pm.corrected
            }
        }else{
            res <- lapply(INDEX, pm.only, x = intensData)
        }
        if(verbose){ 
            cat(" done.\n")
            cat("Computing expression values ...")
        }
        for(k in nr){
            ind <- which(NROW == k)
            temp <- matrix(do.call(rbind, res[ind]), nrow = k)
            ind1 <-  as.vector(sapply(seq_len(n)-1, function(x, ind, m){ ind + x*m }, ind = ind, m = m))
            rob.est[ind1, 1:2] <- robloxbioc(log2(t(temp)), eps = eps, eps.lower = eps.lower, 
                                             eps.upper = eps.upper, steps = steps, mad0 = mad0)
        }
        if(verbose) cat(" done.\n")
        exp.mat <- 2^matrix(rob.est[,1], nrow = m)
        se.mat <- 2^matrix(rob.est[,2], nrow = m)

        dimnames(exp.mat) <- list(ids, sampleNames(x))
        dimnames(se.mat) <- list(ids, sampleNames(x))
        eset <- new("ExpressionSet", phenoData = phenoData(x), 
            experimentData = experimentData(x), exprs = exp.mat, 
            se.exprs = se.mat, annotation = annotation(x))
        return(eset)
    })
