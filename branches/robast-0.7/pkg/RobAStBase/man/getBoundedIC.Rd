\name{getBoundedIC}
\alias{getBoundedIC}

\title{getBoundedIC}
\description{
  Generates a bounded influence curve.
}
\usage{
getBoundedIC(L2Fam, D=trafo(L2Fam@param))
}
\arguments{
  \item{L2Fam}{object of class \code{"L2ParamFamily"}}
  \item{D}{matrix with as many columns as \code{length(L2Fam@param)}}
}
%\details{}
\value{(a bounded) pIC (to matrix \code{D}) given as object of class
 \code{"EuclRandVariable"}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

}
\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
  }
%\note{}
\concept{robust model}
\concept{fixed neighborhood}
\keyword{models}
