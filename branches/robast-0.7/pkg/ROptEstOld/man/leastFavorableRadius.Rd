\name{leastFavorableRadius}
\alias{leastFavorableRadius}
\alias{leastFavorableRadius-methods}
\alias{leastFavorableRadius,L2ParamFamily,UncondNeighborhood,asGRisk-method}

\title{Generic Function for the Computation of Least Favorable Radii}
\description{
  Generic function for the computation of least favorable radii.
}
\usage{
leastFavorableRadius(L2Fam, neighbor, risk, ...)

\S4method{leastFavorableRadius}{L2ParamFamily,UncondNeighborhood,asGRisk}(L2Fam, neighbor, risk, rho, upRad = 1, 
            z.start = NULL, A.start = NULL, upper = 100, maxiter = 100, 
            tol = .Machine$double.eps^0.4, warn = FALSE)
}
\arguments{
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{\dots}{ additional parameters }
  \item{upRad}{ the upper end point of the radius interval to be searched. }
  \item{rho}{ The considered radius interval is: \eqn{[r \rho, r/\rho]}{[r*rho, r/rho]}
    with \eqn{\rho\in(0,1)}{0 < rho < 1}. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
}
%\details{}
\value{
  The least favorable radius and the corresponding inefficiency 
  are computed.
}
\section{Methods}{
\describe{
  \item{L2Fam = "L2ParamFamily", neighbor = "UncondNeighborhood", 
        risk = "asGRisk"}{ computation of the least favorable radius. }
}}
\references{ 
  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Submitted. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{radiusMinimaxIC}}}
\examples{
N <- NormLocationFamily(mean=0, sd=1) 
leastFavorableRadius(L2Fam=N, neighbor=ContNeighborhood(),
                     risk=asMSE(), rho=0.5)
}
\concept{least favorable radius}
\keyword{robust}
