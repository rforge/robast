setMethod("robloxbioc", signature(x = "BeadLevelList"),
    function(x, log = TRUE, imagesPerArray = 1, what = "G", probes = NULL, arrays = NULL,
             eps = NULL, eps.lower = 0, eps.upper = 0.05, steps = 3L, 
             fsCor = TRUE, mad0 = 1e-4){
        BLData <- x
        arraynms <- arrayNames(BLData)
        if(!is.null(arrays) && !is.character(arrays)) arraynms <- arraynms[arrays]
        if(is.character(arrays)) arraynms <- which(arraynms %in% arrays)
        len <- length(arraynms)
        what <- match.arg(what, c("G", "R", "RG", "M", "A", "beta"))
        whatelse <- ""
        if(what == "RG"){
            if(BLData@arrayInfo$channels == "two"){
                what <- "G"
                whatelse <- "R"
            }else{
                stop("Need two-channel data to calculate summary R and G values")
            }
        }
        if(imagesPerArray == 1){
            pr <- getArrayData(BLData, what = "ProbeID", array = arraynms[1])
            sel <- pr != 0
            pr <- pr[sel]
            finten <- getArrayData(BLData, what = what, log = log, array = arraynms[1])[sel]
            nasinf <- !is.finite(finten) | is.na(finten)
            finten <- finten[!nasinf]
        }
        else if(imagesPerArray == 2){
            if(length(arraynms)%%2 != 0) 
                stop("Need an even number of arrays when 'imagesPerArray=2'")
            arrayord <- order(arraynms)
            arraynms <- arraynms[arrayord]
            tmp <- unlist(strsplit(arraynms, "_"))
            chipnums <- tmp[seq(1, length(tmp), by = 3)]
            pos <- tmp[seq(2, length(tmp), by = 3)]
            stripnum <- as.numeric(tmp[seq(3, length(tmp), by = 3)])
            check <- ((chipnums[seq(1, len, by = 2)] == chipnums[seq(2, len, by = 2)]) 
                      & (pos[seq(1, len, by = 2)] == pos[seq(2, len, by = 2)]) 
                      & (stripnum[seq(1, len, by = 2)] == stripnum[seq(2, len, by = 2)] - 1))
            if (sum(check) != length(check)) 
                stop("Missing arrays")
            sel1 <- getArrayData(BLData, what = "ProbeID", array = arraynms[1]) != 0
            sel2 <- getArrayData(BLData, what = "ProbeID", array = arraynms[2]) != 0
            pr <- append(getArrayData(BLData, what = "ProbeID", array = arraynms[1])[sel1], 
                         getArrayData(BLData, what = "ProbeID", array = arraynms[2])[sel2])
            finten <- append(getArrayData(BLData, what = what, log = log, array = arraynms[1])[sel1], 
                             getArrayData(BLData, what = what, log = log, array = arraynms[2])[sel2])
            nasinf <- !is.finite(finten) | is.na(finten)
            finten <- finten[!nasinf]
        }else{
            stop("You can only specify 1 or 2 images per array")
        }
        if(is.null(probes)) probes <- sort(unique(pr))
        probes <- probes[probes > 0 & !is.na(probes)]
        noprobes <- length(probes)
        pr <- pr[!nasinf]
        if (imagesPerArray == 1) {
            G <- GBeadStDev <- GNoBeads <- matrix(0, nrow = noprobes, ncol = len)
            colnames(G) <- colnames(GBeadStDev) <- colnames(GNoBeads) <- arraynms
            if (BLData@arrayInfo$channels == "two" && !is.null(BLData[[arraynms[1]]]$R) && whatelse == "R") 
                R <- RBeadStDev <- RNoBeads <- G
            else R <- NULL
        }
        else if (imagesPerArray == 2) {
            G <- GBeadStDev <- GNoBeads <- matrix(0, nrow = noprobes, ncol = (len/2))
            colnames(G) <- colnames(GBeadStDev) <- colnames(GNoBeads) <- arraynms[seq(1, len, by = 2)]
            if (BLData@arrayInfo$channels == "two" && !is.null(BLData[[arraynms[1]]]$R) && whatelse == "R") 
                R <- RBeadStDev <- RNoBeads <- G
            else R <- NULL
        }
        i <- j <- 1
        while (j <= len) {
            probeIDs <- as.integer(pr)
            start <- 0
            blah <- rmxBeadSummary(x = finten, probeIDs = probeIDs, probes = probes, 
                                   eps = eps, eps.lower = eps.lower, eps.upper = eps.upper, 
                                   steps = steps, fsCor = fsCor, mad0 = mad0)
            G[, i] <- blah$fore
            GBeadStDev[, i] <- blah$sd
            GNoBeads[, i] <- blah$noBeads
            if (BLData@arrayInfo$channels == "two" && !is.null(BLData[[arraynms[i]]]$R) && whatelse == "R") {
                if (imagesPerArray == 1) {
                    finten <- getArrayData(BLData, what = whatelse, log = log, array = arraynms[i])[sel]
                    nasinf <- !is.finite(finten) | is.na(finten)
                    finten <- finten[!nasinf]
                }
                else if (imagesPerArray == 2) {
                    finten <- append(getArrayData(BLData, what = whatelse, log = log, array = arraynms[j])[sel1], 
                                    getArrayData(BLData, what = whatelse, log = log, array = arraynms[j + 1])[sel2])
                    nasinf <- !is.finite(finten) | is.na(finten)
                    finten <- finten[!nasinf]
                }
                blah <- rmxBeadSummary(x = finten, probeIDs = probeIDs, probes = probes, 
                                       eps = eps, eps.lower = eps.lower, eps.upper = eps.upper, 
                                       steps = steps, fsCor = fsCor, mad0 = mad0)
                R[, i] <- blah$fore
                RBeadStDev[, i] <- blah$sd
                RNoBeads[, i] <- blah$noBeads
            }
            j <- j + imagesPerArray
            i <- i + 1
            rm(probeIDs, blah)
            gc()
            if ((imagesPerArray == 1) && (i <= len)) {
                sel <- getArrayData(BLData, what = "ProbeID", array = arraynms[i]) != 0
                pr <- getArrayData(BLData, what = "ProbeID", array = arraynms[i])[sel]
                finten <- getArrayData(BLData, what = what, log = log, array = arraynms[i])[sel]
                nasinf <- !is.finite(finten) | is.na(finten)
                pr <- pr[!nasinf]
                finten <- finten[!nasinf]
            }
            else if ((imagesPerArray == 2) && (j < len)) {
                sel1 <- getArrayData(BLData, what = "ProbeID", array = arraynms[j]) != 0
                sel2 <- getArrayData(BLData, what = "ProbeID", array = arraynms[j + 1]) != 0
                pr <- append(getArrayData(BLData, what = "ProbeID", array = arraynms[j])[sel1], 
                             getArrayData(BLData, what = "ProbeID", array = arraynms[j + 1])[sel2])
                finten <- append(getArrayData(BLData, what = what, log = log, array = arraynms[j])[sel1], 
                                 getArrayData(BLData, what = what, log = log, array = arraynms[j + 1])[sel2])
                nasinf <- !is.finite(finten) | is.na(finten)
                pr <- pr[!nasinf]
                finten <- finten[!nasinf]
            }
        }
        GBeadStDev <- GBeadStDev/sqrt(GNoBeads)
        if(!is.null(R)) RBeadStDev <- RBeadStDev/sqrt(RNoBeads)
        if (whatelse == "R") {
            rownames(G) <- rownames(R) <- rownames(GBeadStDev) <- rownames(RBeadStDev) <- rownames(GNoBeads) <- rownames(RNoBeads) <- probes
            BSData <- new("NChannelSet", R = R, G = G, GBeadStDev = GBeadStDev, 
                          RBeadStDev = RBeadStDev, GNoBeads = GNoBeads, RNoBeads = RNoBeads)
        }
        else {
            BSData <- new("ExpressionSetIllumina")
            assayData(BSData) <- assayDataNew(exprs = G, se.exprs = GBeadStDev, 
                                              NoBeads = GNoBeads, storage.mode = "list")
            rownames(exprs(BSData)) <- rownames(se.exprs(BSData)) <- rownames(NoBeads(BSData)) <- probes
            featureData(BSData) <- new("AnnotatedDataFrame", data = data.frame(ProbeID = probes, row.names = probes))
        }
        if (nrow(pData(BLData)) == len) {
            if (imagesPerArray == 2) 
                BSData@phenoData <- new("AnnotatedDataFrame", data = pData(BLData@phenoData)[arrayord, , drop = FALSE][seq(1, len, by = 2), , drop = FALSE])
            else BSData@phenoData <- BLData@phenoData
        }
        else {
            BSData@phenoData <- new("AnnotatedDataFrame", data = data.frame(sampleName = colnames(G)))
        }
        if (!is.null(pData(BSData)$sampleName)) 
            sampleNames(BSData) <- pData(BSData)$sampleName
        else sampleNames(BSData) <- colnames(G)
        if (whatelse == "R") {
            varMetadata <- data.frame(labelDescription = colnames(BSData@phenoData@data), channel = "_ALL_")
            BSData@phenoData <- new("AnnotatedDataFrame", data = data.frame(BSData@phenoData@data), varMetadata = varMetadata)
        }
        BSData@annotation <- BLData@annotation
        if ("qcScores" %in% names(BLData@arrayInfo)) 
            t <- try(BSData@BeadLevelQC <- BLData@arrayInfo$qcScores, silent = TRUE)
        BSData
    })
## computation of bead summaries via robloxbioc for "matrix"
rmxBeadSummary <- function(x, probeIDs, probes, eps, eps.lower, eps.upper, steps, fsCor, mad0){
    comIDs <- intersect(probeIDs, probes)
    x <- x[probeIDs %in% comIDs]
    probeIDs <- probeIDs[probeIDs %in% comIDs]
    noBeads <- as.vector(table(probeIDs))
    noBeads.uni <- as.integer(names(table(noBeads)))
    probes1 <- comIDs
    len1 <- length(probes1)
    fore1 <- numeric(len1)
    SD1 <- numeric(len1)
    for(i in seq(along = noBeads.uni)){
        index <- noBeads == noBeads.uni[i]
        IDs <- probes1[index]
        if(noBeads.uni[i] == 1){
            fore1[index] <- x[probeIDs %in% IDs]
            SD1[index] <- mad0
        }else{
            temp <- matrix(x[probeIDs %in% IDs], ncol = noBeads.uni[i], byrow = TRUE)
            rmx <- robloxbioc(temp, eps = eps, eps.lower = eps.lower, eps.upper = eps.upper, 
                              steps = steps, fsCor = fsCor, mad0 = mad0)
            fore1[index] <- rmx[,"mean"]
            SD1[index] <- rmx[,"sd"]
        }
    }    
    len <- length(probes)
    fore <- numeric(len)
    SD <- numeric(len)
    noBeads1 <- numeric(len)
    nas <- !(probes %in% comIDs)
    fore[nas] <- NA
    fore[!nas] <- fore1
    SD[nas] <- NA
    SD[!nas] <- SD1
    noBeads1[nas] <- 0
    noBeads1[!nas] <- noBeads

    return(list(fore = fore, sd = SD, noBeads = noBeads1))
}
