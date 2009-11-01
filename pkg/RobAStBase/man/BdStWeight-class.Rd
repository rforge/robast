\name{BdStWeight-class}
\docType{class}
\alias{BdStWeight-class}
\alias{stand,BdStWeight-method}
\alias{stand<-,BdStWeight-method}

\title{Robust Weight classes for bounded, standardized weights}
\description{Classes for bounded, robust, standardized weights.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("BdStWeight", ...)}; 
  to fill slot \code{weight}, you will use the generating functions 
  \code{\link{getweight}} and \code{\link{minbiasweight}}.
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}; inherited from class \code{RobWeight}. }
    \item{\code{weight}}{Object of class \code{"function"} --- the weight function; inherited from class \code{RobWeight}.}
    \item{\code{clip}}{Object of class \code{"numeric"} --- clipping bound(s); inherited from class \code{BoundedWeight}.}
    \item{\code{stand}}{Object of class \code{"matrix"} --- standardization.}
  }
}
\section{Extends}{
Class \code{"RobWeight"}, via class \code{"BoundedWeight"}.
Class \code{"BoundedWeight"}, directly.
}
\section{Methods}{
  \describe{
    \item{stand}{\code{signature(object = "BdStWeight")}:
      accessor function for slot \code{stand}. }

    \item{stand<-}{\code{signature(object = "BdStWeight", value = "matrix")}:
      replacement function for slot \code{stand}. This replacement method 
      should be used with great care, as the slot \code{weight} is not
      simultaneously updated and hence, this may lead to inconsistent
      objects.}
  }
}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}.
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{BoundedWeight-class}}, \code{\link{RobWeight-class}},
         \code{\link{IC}}, \code{\link{InfluenceCurve-class}}}
\examples{
## prototype
new("BdStWeight")
}
\concept{influence curve}
\keyword{classes}
