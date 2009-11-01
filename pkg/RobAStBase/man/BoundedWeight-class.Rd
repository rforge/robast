\name{BoundedWeight-class}
\docType{class}
\alias{BoundedWeight-class}
\alias{clip,BoundedWeight-method}
\alias{clip<-,BoundedWeight-method}

\title{Robust Weight classes for bounded weights}
\description{Classes for bounded, robust weights.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("BoundedWeight", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}; inherited from class \code{RobWeight}. }
    \item{\code{weight}}{Object of class \code{"function"} --- the weight function; inherited from class \code{RobWeight}.}
    \item{\code{clip}}{Object of class \code{"numeric"} --- clipping bound(s).}
  }
}
\section{Extends}{
Class \code{"RobWeight"}, directly.
}
\section{Methods}{
  \describe{
    \item{clip}{\code{signature(x1 = "BoundedWeight")}:
      accessor function for slot \code{clip}. }

    \item{clip<-}{\code{signature(object = "BoundedWeight", value = "numeric")}:
      replacement function for slot \code{clip}. This replacement method 
      should be used with great care, as the slot \code{weight} is not
      simultaneously updated and hence, this may lead to inconsistent
      objects. }
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
\seealso{\code{\link{RobWeight-class}}, \code{\link{IC}}, \code{\link{InfluenceCurve-class}}}
\examples{
## prototype
new("BoundedWeight")
}
\concept{influence curve}
\keyword{classes}
