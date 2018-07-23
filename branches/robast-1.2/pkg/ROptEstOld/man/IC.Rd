\name{IC}
\alias{IC}

\title{Generating function for IC-class}
\description{
  Generates an object of class \code{"IC"}.
}
\usage{
IC(name, Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x}), 
                                        Domain = Reals())), 
   Risks, Infos, CallL2Fam = call("L2ParamFamily"))
}
\arguments{
  \item{name}{ Object of class \code{"character"}. }
  \item{CallL2Fam}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      parametric family. }
  \item{Curve}{ object of class \code{"EuclRandVarList"}. }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link{RiskType-class}}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
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

## The function is currently defined as
IC <- function(name, Curve = EuclRandVarList(RealRandVariable(Map = list(function(x){x})), 
               Domain = Reals()), Risks, Infos, CallL2Fam = call("L2ParamFamily")){
    if(missing(name))
        name <- "square integrable (partial) influence curve"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))
    return(new("IC", name = name, Curve = Curve, Risks = Risks,
               Infos = Infos, CallL2Fam = CallL2Fam))
}
}
\concept{influence curve}
\keyword{robust}
