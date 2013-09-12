\name{CondIC}
\alias{CondIC}

\title{Generating function for CondIC-class}
\description{
  Generates an object of class \code{"CondIC"}.
}
\usage{
CondIC(name, Curve = EuclRandVarList(EuclRandVariable(
       Map = list(function(x){x[1] * x[2]}),
       Domain = EuclideanSpace(dimension = 2),
       Range = Reals())),
       Risks, Infos, CallL2Fam = call("L2RegTypeFamily"))
}
\arguments{
  \item{name}{ character string: name. }
  \item{CallL2Fam}{ object of class \code{"call"}:
      creates an object of \code{"L2RegTypeFamily"}. }
  \item{Curve}{ object of class \code{"EuclRandVariable"}: curve }
  \item{Risks}{ object of class \code{"list"}: 
      list of risks; cf. \code{RiskType-class}. }
  \item{Infos}{ matrix of characters with two columns 
    named \code{method} and \code{message}: additional informations. }
}
%\details{}
\value{Object of class \code{"CondIC"}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}}
\examples{
CondIC()

## The function is currently defined as
function(name, Curve = EuclRandVariable(Map = list(function(x){x[1]*x[2]}),
                             Domain = EuclideanSpace(dimension = 2)),
         Risks, Infos, CallL2Fam = call("L2RegTypeFamily")){
    if(missing(name))
        name <- "Influence curve for a L_2 differentiable regression type family"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))
    return(new("CondIC", name = name, Curve = Curve, Risks = Risks,
               Infos = Infos, CallL2Fam = CallL2Fam))
  }
}
\concept{partial influence curve}
\concept{influence curve}
\keyword{robust}

