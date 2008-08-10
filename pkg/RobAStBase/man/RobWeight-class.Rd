\name{RobWeight-class}
\docType{class}
\alias{RobWeight-class}
\alias{name,RobWeight-method}
\alias{name<-,RobWeight-method}
\alias{weight}
\alias{weight,RobWeight-method}
\alias{weight<-}
\alias{weight<--methods}
\alias{weight<-,RobWeight-method}

\title{Robust Weight classes}
\description{Classes for robust weights.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RobWeight", ...)}. 
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}. }
    \item{\code{weight}:}{Object of class \code{"function"} --- the weight function.}
  }
}
\section{Methods}{
  \describe{
    \item{name}{\code{signature(object = "RobWeight")}: 
      accessor function for slot \code{name}. }

    \item{name<-}{\code{signature(object = "RobWeight")}: 
      replacement function for slot \code{name}. }

    \item{weight}{\code{signature(object = "RobWeight")}: 
      accessor function for slot \code{weight}. }

    \item{weight<-}{\code{signature(object = "RobWeight")}: 
      replacement function for slot \code{weight}. }
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
\seealso{\code{\link{InfluenceCurve-class}}, \code{\link{IC}}}
\examples{
## prototype
new("RobWeight")
}
\concept{influence curve}
\keyword{classes}
