\name{GumbelParameter-class}
\docType{class}
\alias{GumbelParameter-class}
\alias{loc}
\alias{loc,GumbelParameter-method}
\alias{loc<-}
\alias{loc<-,GumbelParameter-method}
\alias{scale,GumbelParameter-method}
\alias{scale<-,GumbelParameter-method}

\title{Paramter of Gumbel distributions}
\description{The class of the parameter of Gumbel distributions.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GumbelParameter", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{loc}}{ real number: location parameter of 
      a Gumbel distribution. }
    \item{\code{scale}}{ positive real number: scale
      parameter of a Gumbel distribution. }
    \item{\code{name}}{ default name is 
      \dQuote{parameter of a Gumbel distribution}. }
  }
}
\section{Extends}{
Class \code{"Parameter"}, directly.\cr
Class \code{"OptionalParameter"}, by class \code{"Parameter"}.
}
\section{Methods}{
  \describe{
    \item{loc}{\code{signature(object = "GumbelParameter")}: access method for
      slot \code{loc}. }
    \item{scale}{\code{signature(x = "GumbelParameter")}: access method for
      slot \code{scale}. }
    \item{loc<-}{\code{signature(object = "GumbelParameter")}: replace method for
      slot \code{loc}. }
    \item{scale<-}{\code{signature(x = "GumbelParameter")}: replace method for
      slot \code{scale}. }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{Gumbel-class}}, \code{\link[distr]{Parameter-class}}}
\examples{
new("GumbelParameter")
}
\concept{Gumbel distribution}
\keyword{distribution}
\concept{parameter}
\concept{S4 parameter class}
\keyword{models}
\concept{generating function}
\concept{scale}
\concept{location}
\concept{location scale model}
