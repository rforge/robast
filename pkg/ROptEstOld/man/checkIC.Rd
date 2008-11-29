\name{checkIC}
\alias{checkIC}

\title{Generic Function for Checking ICs}
\description{
  Generic function for checking centering and Fisher
  consistency of ICs.
}
\usage{
checkIC(IC, L2Fam, ...)
}
\arguments{
  \item{IC}{ object of class \code{"IC"} }
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{\dots}{ additional parameters }
}
\details{ 
  The precisions of the centering and the Fisher consistency
  are computed.
}
\value{The maximum deviation from the IC properties is returned.}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily-class}}, \code{\link{IC-class}}}
\examples{
IC1 <- new("IC")
checkIC(IC1)
}
\concept{influence curve}
\keyword{robust}
