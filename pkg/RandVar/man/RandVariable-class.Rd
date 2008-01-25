\name{RandVariable-class}
\docType{class}
\alias{RandVariable-class}
\alias{Map}
\alias{Domain}
\alias{Range}
\alias{compatibleDomains}
\alias{Map,RandVariable-method}
\alias{Domain,RandVariable-method}
\alias{Range,RandVariable-method}
\alias{Map<-}
\alias{Domain<-}
\alias{Range<-}
\alias{Map<-,RandVariable-method}
\alias{Domain<-,RandVariable-method}
\alias{Range<-,RandVariable-method}
\alias{compatibleDomains,RandVariable,RandVariable-method}
\alias{length,RandVariable-method}
\alias{show,RandVariable-method}

\title{Random variable}
\description{Class of random variables; i.e., measurable maps from 
  \code{Domain} to \code{Range}. The elements contained in 
  the list \code{Map} are functions in one(!) argument named \dQuote{x}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RandVariable", ...)}.
  More frequently they are created via the generating function 
  \code{RandVariable}. 
}
\section{Slots}{
  \describe{
    \item{\code{Map}:}{Object of class \code{"list"}: list of functions. }
    \item{\code{Domain}:}{Object of class \code{"OptionalrSpace"}: 
      domain of the random variable. }
    \item{\code{Range}:}{Object of class \code{"OptionalrSpace"}:
      range of the random variable. }
  }
}
\section{Methods}{
  \describe{
    \item{Map}{\code{signature(object = "RandVariable")}:
      accessor function for the slot \code{Map}.}
    \item{Domain}{\code{signature(object = "RandVariable")}:
      accessor function for the slot \code{Domain}.}
    \item{Range}{\code{signature(object = "RandVariable")}:
      accessor function for the slot \code{Range}.}
    \item{Map<-}{\code{signature(object = "RandVariable")}:
      replacement function for the slot \code{Map}.}
    \item{Domain<-}{\code{signature(object = "RandVariable")}:
      replacement function for the slot \code{Domain}.}
    \item{Range<-}{\code{signature(object = "RandVariable")}:
      replacement function for the slot \code{Range}.}
    \item{compatibleDomains}{\code{signature(e1 = "RandVariable", e2 = "RandVariable")}: 
      test if the domains of two random variables are compatible. }
    \item{length}{\code{signature(object = "RandVariable")}:
      length of the list of functions in slot \code{Map}.}
    \item{show}{\code{signature(object = "RandVariable")}}
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{RandVariable}}, \code{\link{EuclRandVariable-class}},
  \code{\link{EuclRandMatrix-class}}, \code{\link{EuclRandVarList-class}}}
\examples{
(R1 <- new("RandVariable"))
Map(R1)
Domain(R1)
Range(R1)
Map(R1) <- list(function(x){ceiling(x)}, function(x){floor(x)})
Domain(R1) <- Reals()
Range(R1) <- Naturals()
R1
Map(R1)
length(R1)

R2 <- R1
Domain(R2) <- Naturals()
compatibleDomains(R1, R2)
Domain(R2) <- NULL
compatibleDomains(R1, R2)
Domain(R2) <- EuclideanSpace(dimension = 1)
compatibleDomains(R1, R2)
Domain(R2) <- EuclideanSpace(dimension = 2)
compatibleDomains(R1, R2)
}
\concept{random variable}
\keyword{classes}
