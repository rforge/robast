\name{InfluenceCurve}
\alias{InfluenceCurve}

\title{Generating function for InfluenceCurve-class}
\description{
  Generates an object of class \code{"InfluenceCurve"}.
}
\usage{
InfluenceCurve(name, Curve = EuclRandVarList(EuclRandVariable(Domain = Reals())), 
            Risks, Infos)
}
\arguments{
  \item{name}{ character string: name of the influence curve }
  \item{Curve}{ object of class \code{"EuclRandVarList"} }
  \item{Risks}{ list of risks  }
  \item{Infos}{ matrix of characters with two columns
    named \code{method} and \code{message}: additional informations }
}
%\details{}
\value{Object of class \code{"InfluenceCurve"}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.
  
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve-class}}}
\examples{
InfluenceCurve()

## The function is currently defined as
InfluenceCurve <- function(name, Curve = EuclRandVarList(EuclRandVariable(Domain = Reals())), 
                           Risks, Infos){
    if(missing(name))
        name <- "influence curve"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))
    
    return(new("InfluenceCurve", name = name, Curve = Curve, 
               Risks = Risks, Infos = Infos))
}
}
\concept{influence curve}
\keyword{robust}
