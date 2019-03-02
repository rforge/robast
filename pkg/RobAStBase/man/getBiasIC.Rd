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
             biastype = symmetricBias(), normtype = NormType(),
             tol = .Machine$double.eps^0.25, numbeval = 1e5, withCheck = TRUE, ...)
}
\arguments{
  \item{IC}{ object of class \code{"InfluenceCurve"} }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"}. }
  \item{biastype}{object of class \code{"BiasType"}}
  \item{normtype}{object of class \code{"NormType"}}
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{numbeval}{number of evalation points.}
  \item{withCheck}{logical: should a call to \code{checkIC} be done to
                   check accuracy (defaults to \code{TRUE}).}
  \item{\dots}{ additional parameters to be passed to expectation \code{E} }
}
%\details{}
\value{The bias of the IC is computed.}
\section{Methods}{
\describe{
  \item{IC = "IC", neighbor = "UncondNeighborhood"}{
    determines the as. bias by random evaluation of the IC;
    this random evaluation is done by the internal S4-method
    \code{.evalBiasIC}; this latter dispatches according to
    the signature \code{IC, neighbor, biastype}.\cr
    For signature \code{IC="IC", neighbor = "ContNeighborhood",
        biastype = "BiasType"}, also an argument \code{normtype}
    is used to be able to use self- or information standardizing
    norms; besides this the signatures
    \code{IC="IC", neighbor = "TotalVarNeighborhood",
        biastype = "BiasType"},
        \code{IC="IC", neighbor = "ContNeighborhood",
        biastype = "onesidedBias"}, and
    \code{IC="IC", neighbor = "ContNeighborhood",
        biastype = "asymmetricBias"} are implemented.
    }
}}
\references{
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.

  Ruckdeschel, P. and Kohl, M. (2005) Computation of the Finite Sample Bias
  of M-estimators on Neighborhoods.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\note{This generic function is still under construction.}
\seealso{\code{\link{getRiskIC-methods}}, \code{\link{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
