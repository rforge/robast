\name{rgsOptIC.ALs}
\alias{rgsOptIC.ALs}

\title{Computation of the optimally robust IC for ALs estimators}
\description{
  The function \code{rgsOptIC.ALs} computes the optimally robust IC
  for ALs estimators in case of linear regression with unknown
  scale and (convex) contamination neighborhoods where the
  regressor is random; confer Subsection 7.3.1 of Kohl (2005).
}
\usage{
rgsOptIC.ALs(r, K, A.rg.start, b.rg.Up = 1000, delta = 1e-06, 
             itmax = 50, check = FALSE)
}
\arguments{
  \item{r}{ non-negative real: neighborhood radius. }
  \item{K}{ object of class \code{"Distribution"}. }
  \item{A.rg.start}{ positive definite and symmetric matrix:
    starting value for the standardizing matrix of the
    regression part. }
  \item{b.rg.Up}{ positive real: the upper end point of the
    interval to be searched for b.rg. }
  \item{delta}{ the desired accuracy (convergence tolerance). }
  \item{itmax}{ the maximum number of iterations. }
  \item{check}{ logical. Should constraints be checked. }
}
\details{
  If \code{A.rg.start} is missing, the inverse of the 
  second moment matrix of \code{K} is used.
}
\value{Object of class \code{"ContIC"}}
\references{ 
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{ContIC-class}}
\examples{
K <- DiscreteDistribution(1:5) # = Unif({1,2,3,4,5})
IC1 <- rgsOptIC.ALs(r = 0.1, K = K)
checkIC(IC1)
Risks(IC1)
Infos(IC1)
}
\concept{regression and scale}
\concept{influence curve}
\keyword{robust}
