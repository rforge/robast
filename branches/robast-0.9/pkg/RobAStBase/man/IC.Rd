\name{IC}
\alias{IC}

\title{Generating function for IC-class}
\description{
  Generates an object of class \code{"IC"}.
}
\usage{
IC(name, Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                        Domain = Reals())), 
   Risks, Infos, CallL2Fam = call("L2ParamFamily"), modifyIC = NULL)
}
\arguments{
  \item{name}{ Object of class \code{"character"}; the name of the IC. }
  \item{CallL2Fam}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      parametric family. }
  \item{Curve}{ object of class \code{"EuclRandVarList"}. }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link[distrMod]{RiskType-class}}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
  \item{modifyIC}{ Object of class \code{"OptionalFunction"}:
    function of two arguments, which are an L2 parametric family 
    and an optional influence curve. Returns an object of 
    class \code{"IC"} at the parameter value of the L2 parametric
    family. This function is mainly used for internal
    computations! }
}
%\details{}
\value{Object of class \code{"IC"}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{IC-class}}}
\examples{
IC1 <- IC()
plot(IC1)
}
\concept{influence curve}
\keyword{robust}
