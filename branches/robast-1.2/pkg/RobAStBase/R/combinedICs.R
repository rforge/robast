################################################################################
if(FALSE){
################################################################################
## 20180809: reverted changes from rev 1110
################################################################################

combineOrthPICs <- function(pIC1, pIC2, combinedName = "combined IC", dim){
  ## adds to complementary pICs to give one IC
  ## the orthogonality is not checked here

    IC <- new(".fastIC")
    IC@name <- combinedName
    pICC1 <- as(diag(dim)%*%pIC1@Curve,"EuclRandVariable")
    pICC2 <- as(diag(dim)%*%pIC2@Curve,"EuclRandVariable")
    IC@Curve <- EuclRandVarList(pICC1+pICC2)
    IC@Risks <- pIC1@Risks
    if(length(pIC2@Risks)) addRisk(IC) <- pIC2@Risks
    IC@Infos <- pIC1@Infos
    if(nrow(pIC2@Infos)) addInfo(IC) <- pIC2@Infos
    IC@CallL2Fam <- pIC1@CallL2Fam
    .modifyIC.0 <- function(L2Fam, IC, withMakeIC = FALSE){
       pic1 <- pic1@modifyIC(L2Fam, pIC1, withMakeIC)
       pic2 <- pic2@modifyIC(L2Fam, pIC2, withMakeIC)
       IC1 <- combineOrthPICs(pic1, pic2,combinedName)
       return(IC1)
    }
    .modifyIC.1 <- function(L2Fam, IC, withMakeIC = FALSE){
       IC1 <- .modifyIC.0(L2Fam, IC, withMakeIC)
       IC1@modifyIC <- .modifyIC.1
       return(IC1)
    }

    IC@modifyIC <- .modifyIC.1
    IC@.fastFct <- function(x){pIC1@.fastFct(x)+pIC2@.fastFct(x)}
    return(IC)
}


.fastIC <- function(name ="", Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}),
               Domain = Reals())), Risks, Infos, CallL2Fam = call("L2ParamFamily"),
               modifyIC = NULL, .fastFct = NULL){
fastIC <- new(".fastIC")
if(missing(Infos)) Infos <- fastIC@Infos
if(missing(Risks)) Risks <- fastIC@Risks
IC.0 <- IC(name, Curve, Risks, Infos, CallL2Fam, modifyIC)
slotNms <- slotNames(class(IC.0))
for(sN in slotNms) slot(fastIC, sN) <- slot(IC.0,sN)
if(is.null(.fastFct)||missing(.fastFct)){
   ICM <- IC.0@Curve[[1]]@Map
    .fastFct <- function(x){
       if(is.null(dim(x)))
          sapply(x, function(u) sapply(ICM, function(s)s(u)))
       else
          apply(x, 1,function(u) sapply(ICM, function(s)s(u)))
    }
}
fastIC@.fastFct <- .fastFct
return(fastIC)
}
################################################################################
## end if(FALSE)
################################################################################
}