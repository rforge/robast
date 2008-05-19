\name{InfRobRegTypeModel-class}
\docType{class}
\alias{InfRobRegTypeModel-class}
\alias{neighbor<-,InfRobRegTypeModel-method}
\alias{show,InfRobRegTypeModel-method}

\title{Robust regression-type model with infinitesimal neighborhood}
\description{Class of robust regression-type models with infinitesimal 
  (conditional or unconditional) neighborhoods; i.e., the neighborhood 
  is shrinking at a rate of \eqn{\sqrt{n}}{sqrt(n)}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("InfRobRegTypeModel", ...)}.
  More frequently they are created via the generating function 
  \code{InfRobRegTypeModel}. 
}
\section{Slots}{
  \describe{
    \item{\code{center}:}{Object of class \code{"L2RegTypeFamily"}. }
    \item{\code{neighbor}:}{Object of class \code{"Neighborhood"}. }
  }
}
\section{Extends}{
Class \code{"RobModel"}, directly.
}
\section{Methods}{
  \describe{
    \item{neighbor<-}{\code{signature(object = "InfRobRegTypeModel")}: 
      replacement function for slot \code{neighbor}.}
    \item{show}{\code{signature(object = "InfRobRegTypeModel")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2RegTypeFamily-class}}, \code{Neighborhood-class},
  \code{\link{InfRobRegTypeModel}}}
\examples{
new("InfRobRegTypeModel")
}
\concept{robust model}
\concept{infinitesimal neighborhood}
\concept{shrinking neighborhood}
\keyword{classes}
\keyword{models}
