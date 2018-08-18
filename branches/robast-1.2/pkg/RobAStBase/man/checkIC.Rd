\name{checkIC}
\alias{checkIC}
\alias{checkIC,IC,missing-method}
\alias{checkIC,IC,L2ParamFamily-method}

\title{Generic Function for Checking ICs}
\description{
  Generic function for checking centering and Fisher
  consistency of ICs.
}
\usage{
checkIC(IC, L2Fam, ...)
\S4method{checkIC}{IC,missing}(IC, out = TRUE, ..., diagnostic = FALSE)
\S4method{checkIC}{IC,L2ParamFamily}(IC, L2Fam, out = TRUE,..., diagnostic = FALSE)
}
\arguments{
  \item{IC}{ object of class \code{"IC"} }
  \item{L2Fam}{ L2-differentiable family of probability measures. }
  \item{out}{ logical: Should the values of the checks be printed out?}
  \item{\dots}{ additional parameters }
  \item{diagnostic}{ logical; if \code{TRUE} and \code{out==TRUE},
    diagnostic information on the integration is printed; independent
    of \code{out}, if \code{diagnostic==TRUE}, this information is
    returned as attribute \code{diagnostic} of the return value. . }
}
\details{ 
  The precisions of the centering and the Fisher consistency
  are computed.

  Diagnostics on the involved integrations are available if argument
  \code{diagnostic} is \code{TRUE}. Then there is attribute \code{diagnostic}
  attached to the return value, which may be inspected
  and accessed through \code{\link[distrEx:distrExIntegrate]{showDiagnostic}} and
  \code{\link[distrEx:distrExIntegrate]{getDiagnostic}}.
}
\value{The maximum deviation from the IC properties is returned.}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}}
\examples{
IC1 <- new("IC")
checkIC(IC1)
}
\concept{influence curve}
\keyword{robust}
