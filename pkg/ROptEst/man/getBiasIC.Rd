\name{getBiasIC}
\alias{getBiasIC}
\alias{getBiasIC-methods}
\alias{getBiasIC,IC,UncondNeighborhood-method}

\title{Generic function for the computation of the asymptotic bias for an IC}
\description{
  Generic function for the computation of the asymptotic bias for an IC.
}
\usage{
getBiasIC(IC, neighbor, ...)

\S4method{getBiasIC}{IC,UncondNeighborhood}(IC, neighbor, L2Fam, 
                     biastype = symmetricBias(), 
                     tol = .Machine$double.eps^0.25)

}
\arguments{
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{\dots}{ additional parameters }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"} or missing. }
  \item{biastype}{ object of class \code{"BiasType"}. }
  \item{tol}{ the desired accuracy (convergence tolerance).}  
}
\details{To make sure that the results are valid, it is recommended
  to include an additional check of the IC properties of \code{IC} 
  using \code{checkIC}.}
\value{The asymptotic bias of an IC is computed.}
\section{Methods}{
\describe{

  \item{IC = "IC", neighbor = "UncondNeighborhood":}{
   asymptotic bias of \code{IC} in case of unconditional neighborhoods. }

}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},
        Peter Ruckdeschel \email{Peter.Ruckdeschel@uni-bayreuth.de}}

\note{This generic function is still under construction.}
\seealso{\code{\link{getRiskIC-methods}}, \code{\link[RobAStBase]{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{}
