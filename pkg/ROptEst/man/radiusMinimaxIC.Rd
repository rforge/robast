\name{radiusMinimaxIC}
\alias{radiusMinimaxIC}
\alias{radiusMinimaxIC-methods}
\alias{radiusMinimaxIC,L2ParamFamily,UncondNeighborhood,asGRisk-method}

\title{Generic function for the computation of the radius minimax IC}
\description{
  Generic function for the computation of the radius minimax IC.
}
\usage{
radiusMinimaxIC(L2Fam, neighbor, risk, ...)

\S4method{radiusMinimaxIC}{L2ParamFamily,UncondNeighborhood,asGRisk}(
        L2Fam, neighbor, risk, loRad, upRad, z.start = NULL, A.start = NULL, 
        upper = 1e5, maxiter = 50, tol = .Machine$double.eps^0.4, 
        warn = FALSE, verbose = FALSE)
}
\arguments{
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{\dots}{ additional parameters. }
  \item{loRad}{ the lower end point of the interval to be searched. }
  \item{upRad}{ the upper end point of the interval to be searched. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed }
}
\details{
In case the neighborhood radius is unknown, Rieder et al. (2001, 2008)
and Kohl (2005) show that there is nevertheless a way to compute an
optimally robust IC - the so-called radius-minimax IC - which is 
optimal for some radius interval.
}
\value{The radius minimax IC is computed.}
\section{Methods}{
\describe{
  \item{L2Fam = "L2ParamFamily", neighbor = "UncondNeighborhood", risk = "asGRisk":}{ 
    computation of the radius minimax IC for an L2 differentiable parametric family. }
}}
\references{ 
  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications, \emph{17}(1) 13-40.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{radiusMinimaxIC}}}
\examples{
N <- NormLocationFamily(mean=0, sd=1) 
radIC <- radiusMinimaxIC(L2Fam=N, neighbor=ContNeighborhood(), 
                         risk=asMSE(), loRad=0.1, upRad=0.5)
checkIC(radIC)
}
\concept{radius minimax influence curve}
\concept{influence curve}
\keyword{robust}
