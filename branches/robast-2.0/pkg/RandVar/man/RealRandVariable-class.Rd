\name{RealRandVariable-class}
\docType{class}
\alias{RealRandVariable-class}
\alias{Range<-,RealRandVariable-method}

\title{Real random variable}
\description{Class of real random variables.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RealRandVariable", ...)}.
  More frequently they are created via the generating function 
  \code{EuclRandVariable}. 
}
\section{Slots}{
  \describe{
    \item{\code{Map}:}{Object of class \code{"list"}: list of functions. }
    \item{\code{Domain}:}{Object of class \code{"OptionalrSpace"}: 
      domain of the random variable. }
    \item{\code{Range}:}{Object of class \code{"Reals"}:
      range of the random variable. }
  }
}
\section{Extends}{
Class \code{"EuclRandVariable"}, directly.\cr
Class \code{"RandVariable"}, by class \code{"EuclRandVariable"}.
}
\section{Methods}{
  \describe{
    \item{Range<-}{\code{signature(object = "EuclRandVariable")}:
      replacement function for the slot \code{Range}. }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{EuclRandVariable-class}}}
\examples{
new("RealRandVariable", Map=list(function(x){x}), Range = Reals())
}
\concept{random variable}
\keyword{classes}
