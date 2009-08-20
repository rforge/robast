\name{HampelWeight-class}
\docType{class}
\alias{HampelWeight-class}
\alias{cent,HampelWeight-method}
\alias{cent<-,HampelWeight-method}

\title{Robust Weight classes for weights of Hampel type}
\description{Classes for weights of Hampel type.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("HampelWeight", ...)};
  to fill slot \code{weight}, you will use the generating functions 
  \code{\link{getweight}} and \code{\link{minbiasweight}}.
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}; inherited from class \code{RobWeight}. }
    \item{\code{weight}}{Object of class \code{"function"} --- the weight function; inherited from class \code{RobWeight}.}
    \item{\code{clip}}{Object of class \code{"numeric"} --- clipping bound(s); inherited from class \code{BoundedWeight}.}
    \item{\code{stand}}{Object of class \code{"matrix"} --- standardization; inherited from class \code{BdStWeight}.}
    \item{\code{cent}}{Object of class \code{"numeric"} --- centering.}
  }
}
\section{Extends}{
Class \code{"RobWeight"}, via class \code{"BoundedWeight"}.
Class \code{"BoundedWeight"}, via class \code{"BdStWeight"}.
Class \code{"BdStWeight"}, directly.
}
\section{Methods}{
  \describe{
    \item{cent}{\code{signature(object = "HampelWeight")}:
      accessor function for slot \code{cent}. }

    \item{cent<-}{\code{signature(object = "HampelWeight", value = "matrix")}:
      replacement function for slot \code{cent}. This replacement method 
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
\seealso{\code{\link{BdStWeight-class}},
         \code{\link{BoundedWeight-class}}, \code{\link{RobWeight-class}},
         \code{\link{IC}}, \code{\link{InfluenceCurve-class}}}
\examples{
## prototype
new("HampelWeight")
}
\concept{influence curve}
\keyword{classes}
