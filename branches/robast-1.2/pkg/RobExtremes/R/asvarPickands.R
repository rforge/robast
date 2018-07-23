asvarPickands <- function(model, alpha=2){

    isGP <- is(model,"GParetoFamily")
    if(!(isGP|is(model,"GEVFamily")))
         stop("Pickands estimator only available for GPD and GEVD.")

  scshn <- scaleshapename(model)
#  par0 <- main(model@param)[scshn]
#  beta <- par0[1]; xi <- par0[2]

  if(isGP){
    al1 <- 1-1/alpha
    al2 <- 1-1/alpha^2
  }else{
    al1 <- exp(-1/alpha)
    al2 <- exp(-1/alpha^2)
  }

  M2 <- q.l(model)(al1)
  M4 <- q.l(model)(al2)

  xi <- log((M4-M2)/M2)/log(alpha)
  qu <- 1/(alpha^xi-1)
  beta <- xi * M2 * qu
  # d/dMi xi = h11, h12
  # d/dMi beta = (d/dMi M2) * xi* qu + (d/dMi xi) * M2 * qu +(d/dMi qu)* M2 * xi
  # d/dMi M2 = 1 for i=2, = 0 for i=4
  # d/dMi qu = d/dxi qu * d/dMi xi
  # d/dxi qu = - qu^2 * alpha^xi * log(alpha)
  # => d/dMi beta = (i==2)*xi*qu +
  #        (M2 * qu - beta * qu * alpha^xi * log(alpha)) * (h11,h12)
  # dqu =  M2 * qu + beta (d/dxi qu /qu)
  dqu <-  M2 * qu - beta * qu * alpha^xi * log(alpha)
  h11 <- -M4/(M2*(M4-M2))/log(alpha)
  h12 <- 1/(M4-M2)/log(alpha)

  h21 <- h11*dqu + xi*qu
  h22 <- h12*dqu

  ### corresponding terms for original definition for beta, i.e.
  ##   beta <- xi * M2^2/(M4-2*M2)
  #t1 <- 2*M2*(M4-M2)/(M4-2*M2)^2
  #t2 <- -M2^2/(M4-2*M2)^2
  #h21 <- h11*M2^2/(M4-2*M2) + t1*log((M4-M2)/M2)/log(alpha)
  #h22 <- h12*M2^2/(M4-2*M2) + t2*


  C <- matrix(c(h21,h22,h11,h12),2,2)

#  f1 <- (1-al1)^(1+xi)/beta
#  f2 <- (1-al2)^(1+xi)/beta
#  M <- matrix(c(al1-1,al2-1,al1,al2-1,al1,al2),ncol=3)
#  Werte <- t(C) %*% diag(1/c(f1,f2)) %*% M
#  GES <- max(colSums(Werte^2)^.5)
#  GES

  if(isGP){
  s11 <- al1*(1-al1)^(-1-2*xi)
  s12 <- al1*(1-al1)^(-1-xi)*(1-al2)^(-xi)
  s21 <- s12
  s22 <- al2*(1-al2)^(-1-2*xi)
  }else{
  s11 <- al1^(-1)*(1-al1)*(-log(al1))^(-2-2*xi)
  s12 <- al2^(-1)*(1-al2)*(log(al1)*log(al2))^(-1-1*xi)
  s21 <- s12
  s22 <- al2^(-1)*(1-al2)*(-log(al2))^(-2-2*xi)
  }
  S <- beta^2*matrix(c(s11,s12,s21,s22),2,2)

  ASV_Pick <- t(C) %*% S %*% (C)
  ASV_Pick <- PosSemDefSymmMatrix(ASV_Pick)
  dimnames(ASV_Pick) <- list(scshn,scshn)
  return(ASV_Pick)
}

asvarQBCC <- function(model, p1 = 1/3, p2= 2/3){

   if(!(is(model,"WeibullFamily")))
         stop("QuantileBCC estimator only available for Weibull.")

  scshn <- scaleshapename(model)

 if(p1>=p2) {p<-p1; p1 <- p2; p2 <- p}

  qm <- q.l(model)
  Q1 <- qm(p1)
  Q2 <- qm(p2)
  l1 <- -log(p1); l2 <- -log(p2)

  lq <- 1/(log(l2)-log(l1))
  xi <- (log(Q2)-log(Q1))*lq
  beta <- Q2*l2^(-1/xi)

  dqu <-  beta * log(l2) /xi^2
  h11 <- -lq/Q1
  h12 <- lq/Q2

  h21 <- h11*dqu
  h22 <- h12*dqu + l2^(-1/xi)

  C <- matrix(c(h21,h22,h11,h12),2,2)
  dm <- d(model)
  s11 <- p1*(1-p1)/dm(Q1)^2
  s12 <- p1*(1-p2)/dm(Q1)/dm(Q2)
  s21 <- s12
  s22 <- p2*(1-p2)/dm(Q2)^2

  S <- matrix(c(s11,s12,s21,s22),2,2)

  ASV <- t(C) %*% S %*% (C)
  ASV <- PosSemDefSymmMatrix(ASV)
  dimnames(ASV) <- list(scshn,scshn)
  return(ASV)
}











