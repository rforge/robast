\name{FixRobRegTypeModel-class}
\docType{class}
\alias{FixRobRegTypeModel-class}
\alias{neighbor<-,FixRobRegTypeModel-method}
\alias{show,FixRobRegTypeModel-method}

\title{Robust regression-type model with fixed neighborhood}
\description{Class of robust regression-type models with fixed 
  (conditional or unconditional) neighborhoods.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("FixRobRegTypeModel", ...)}.
  More frequently they are created via the generating function 
  \code{FixRobRegTypeModel}. 
}
\section{Slots}{
  \describe{
    \item{\code{center}:}{Object of class \code{"RegTypeFamily"}. }
    \item{\code{neighbor}:}{Object of class \code{"Neighborhood"}. }
  }
}
\section{Extends}{
Class \code{"RobModel"}, directly.
}
\section{Methods}{
  \describe{
    \item{neighbor<-}{\code{signature(object = "FixRobRegTypeModel")}
      replacement function for slot \code{neighbor}.}
    \item{show}{\code{signature(object = "FixRobRegTypeModel")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{RegTypeFamily-class}}, \code{\link[ROptEstOld]{Neighborhood-class}},
  \code{\link{FixRobRegTypeModel}}}
\examples{
new("FixRobRegTypeModel")
}
\concept{robust model}
\concept{fixed neighborhood}
\keyword{classes}
\keyword{models}
\keyword{robust}
