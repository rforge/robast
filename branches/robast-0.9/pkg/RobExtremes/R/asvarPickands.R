asvarPickands <- function(model, alpha=2){

    isGP <- is(ParamFamily,"GParetoFamily")
    if(!(isGP|is(ParamFamily,"GEVFamily")))
         stop("Pickands estimator only available for GPD and GEVD.")

  scshn <- scaleshapename(model)
  par0 <- main(model@param)[scshn]
  beta <- par0[1]; xi <- par0[2]

  if(isGP){
    al1 <- 1-1/alpha
    al2 <- 1-1/alpha^2
  }else{
    al1 <- exp(-1/alpha)
    al2 <- exp(-1/alpha^2)
  }
  M2 <- q(model)(al1)
  M4 <- q(model)(al2)

  h11 <- -M4/(M2*(M4-M2))/log(alpha)
  h12 <- 1/(M4-M2)/log(alpha)
  t1 <- 2*M2*(M4-M2)/(M4-2*M2)^2
  t2 <- -M2^2/(M4-2*M2)^2
  h21 <- h11*M2^2/(M4-2*M2) + t1*log((M4-M2)/M2)/log(alpha)
  h22 <- h12*M2^2/(M4-2*M2) + t2*log((M4-M2)/M2)/log(alpha)


  C <- matrix(c(h21,h22,h11,h12),2,2)

#  f1 <- (1-al1)^(1+xi)/beta
#  f2 <- (1-al2)^(1+xi)/beta
#  M <- matrix(c(al1-1,al2-1,al1,al2-1,al1,al2),ncol=3)
#  Werte <- t(C) %*% diag(1/c(f1,f2)) %*% M
#  GES <- max(colSums(Werte^2)^.5)
#  GES

  s11 <- al1*(1-al1)^(-1-2*xi)
  s12 <- al1*(1-al1)^(-1-xi)*(1-al2)^(-xi)
  s21 <- s12
  s22 <- al2*(1-al2)^(-1-2*xi)

  S <- beta^2*matrix(c(s11,s12,s21,s22),2,2)

  ASV_Pick <- t(C) %*% S %*% (C)
  ASV_Pick <- PosSemDefSymmMatrix(ASV_Pick)
  dimnames(ASV_Pick) <- list(scshn,scshn)
  return(ASV_Pick)
}












