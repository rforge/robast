\name{getBoundedIC}
\alias{getBoundedIC}

\title{getBoundedIC}
\description{
  Generates a bounded influence curve.
}
\usage{
getBoundedIC(L2Fam, D=trafo(L2Fam@param), ..., diagnostic = FALSE)
}
\arguments{
  \item{L2Fam}{object of class \code{"L2ParamFamily"}}
  \item{D}{matrix with as many columns as \code{length(L2Fam@param)}}
  \item{...}{further arguments to be passed to \code{E}}
  \item{diagnostic}{ logical; if \code{TRUE}, the return value obtains
    an attribute \code{"diagnostic"} with diagnostic information on the
    integration. }
}
%\details{}
\value{(a bounded) pIC (to matrix \code{D}) given as object of class
 \code{"EuclRandVariable"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }
%\note{}
\concept{robust model}
\concept{fixed neighborhood}
\keyword{models}
