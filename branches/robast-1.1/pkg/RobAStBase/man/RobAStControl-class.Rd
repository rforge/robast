\name{RobAStControl-class}
\docType{class}
\alias{RobAStControl-class}
\alias{name,RobAStControl-method}
\alias{name<-,RobAStControl-method}

\title{Control classes in package RobAStBase}
\description{Control classes in package \pkg{RobAStBase}.}
\section{Objects from the Class}{
  This class is virtual; that is no objects may be created.
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}:
      name of the control object. }
  }
}
\section{Methods}{
  \describe{
    \item{name}{\code{signature(object = "RobAStControl")}: 
      accessor function for slot \code{name}. }

    \item{name<-}{\code{signature(object = "RobAStControl", value = "character")}: 
      replacement function for slot \code{name}. }
  }
}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
%\seealso{}
%\examples{}
\concept{influence curve}
\keyword{classes}
