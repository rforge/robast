gam3 <- function(x) gamma(x)*(psigamma(x,3)+4*psigamma(x,2)*digamma(x)+3*trigamma(x)^2+
                              6*digamma(x)^2*trigamma(x)+digamma(x)^4)
gam2 <- function(x) gamma(x)*(psigamma(x,2)+3*trigamma(x)*digamma(x)+digamma(x)^3)
gam1 <- function(x) gamma(x)*(trigamma(x)+digamma(x)^2)
gam0 <- function(x) gamma(x)*digamma(x)


I11 <- 1
I12 <- -gam0(3)+2*gam0(2)-gam0(1)

I22 <- gam1(1)-2*gam1(2)+gam1(3)+2*gam0(1)-2*gam0(2)+1

I13 <- -gam0(2)+gam1(3)/2-gam1(2)/2

I23 <- gam2(3)/2-gam2(2)/2+gam1(2)-gam1(1)

I33 <- (gam3(3)-2*gam3(2)+gam3(1))/4+gam2(1)-gam2(2)+gam1(1)

(Ima <-matrix(c(I11,I12,I13,I12,I22,I23,I13,I23,I33),3,3))

check <- function(a=3,b=0.8){
  ga1 <- gamma(a+b)
  ga2 <- c(gamma(a),
           gamma(a)+gam0(a)*b,
           gamma(a)+gam0(a)*b+gam1(a)*b^2/2,
           gamma(a)+gam0(a)*b+gam1(a)*b^2/2+gam2(a)*b^3/6,
           gamma(a)+gam0(a)*b+gam1(a)*b^2/2+gam2(a)*b^3/6+gam3(a)*b^4/24
           )
  print(c(ga1,ga2,ga1-ga2))
}
