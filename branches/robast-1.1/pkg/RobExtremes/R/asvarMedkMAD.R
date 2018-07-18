asvarMedkMAD <- function( model, k=1){
  if(! is(model, "L2ScaleShapeUnion"))
     stop("This function only works for Scale-Shape models")

  scshn <- scaleshapename(model)
  par0 <- main(model@param)[scshn]
  beta <- par0[1]; xi <- par0[2]

  M <- kMAD(model@distribution, k=k)
  m <- q.l(model)(.5)

  x1.0 <- m - M
  x2.0 <- m + k * M

  ## joint Variance of median and kMAD, see Serfling Mazumder
  dmm <- d(model)(x1.0)
  dmp <- d(model)(x2.0)
  dm <-  d(model)(m)
  alpha <- p(model)(m-M)+p(model)(m+k*M)
  betA <- dmm-dmp
  ceta <- dmm+k*dmp
  eta <- betA^2 + 4*(1-alpha)*betA*dm

  g22 <- 1/4/dm^2
  g12 <- 1/4/dm/ceta*(1-4*p(model)(m-M)+betA/dm)
  g11 <- 1/4/ceta^2*(1+eta/dm^2)

  V <- matrix(c(g11,g12,g12,g22),2,2)

  if(is(model,"GParetoFamily")){
      uf <- function(x) 1+xi*x/beta
      u1.0 <- uf(x1.0)
      u2.0 <- uf(x2.0)
      um.0 <- uf(m)

      dfc <- function(u1,u2,s1,s2,fct) s2*fct(u2)-s1*fct(u1)

      uxi <- function(u) u^(-1/xi-1)

      Gxi <- function(u) uxi(u)*(u*(1-log(u))-1)/xi^2
      Gbet <- function(u) uxi(u)*(1-u)/xi/beta
      Gx <- function(u) uxi(u)/beta

      dG1_xi  <- dfc(u1=u1.0,u2=u2.0,s1= 1,s2=1,fct=Gxi)
      dG1_beta <- dfc(u1=u1.0,u2=u2.0,s1= 1,s2=1,fct=Gbet)
      dG1_M <-   dfc(u1=u1.0,u2=u2.0,s1=-1,s2=k,fct=Gx)
      dG1_m <-   dfc(u1=u1.0,u2=u2.0,s1= 1,s2=1,fct=Gx)

      dG2_xi  <- Gxi(um.0)
      dG2_beta <- Gbet(um.0)
      dG2_M <- 0
      dG2_m <- Gx(um.0)

      D1 <- matrix(c(dG1_beta,dG2_beta,dG1_xi,dG2_xi),2,2)
      D2 <- matrix(c(dG1_M,dG2_M,dG1_m,dG2_m),2,2)

      D <- -solve(D1)%*%D2
  }else{
   psi_med <- function(x) (0.5-(x<=m))/dm
   psi_kMad <- function(x){
       cp <- k*dmp+dmm
       cm = dmp-dmm
       return((0.5-((x<=m+k*M)&(x>=m-M)))/cp + cm/cp*((x<=m)-0.5)/dm)
   }

   L2d <- model@L2deriv[[1]]
   L_xi.f = function(x) evalRandVar(L2d,x)[2,]
   L_beta.f = function(x) evalRandVar(L2d,x)[1,]

   E11 <- E(distribution(model),fun=function(x) psi_kMad(x) * L_beta.f(x))
   E12 <- E(distribution(model),fun=function(x) psi_kMad(x) * L_xi.f(x))
   E21 <- E(distribution(model),fun=function(x) psi_med(x) * L_beta.f(x))
   E22 <- E(distribution(model),fun=function(x) psi_med(x) * L_xi.f(x))
   D <- solve(matrix(c(E11,E21,E12,E22),2,2))
  }

  ASV_Med <- PosSemDefSymmMatrix(D %*% V %*% t(D))
  dimnames(ASV_Med) <- list(scshn,scshn)
  return(ASV_Med)
}

