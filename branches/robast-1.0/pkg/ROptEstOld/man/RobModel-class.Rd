\name{RobModel-class}
\docType{class}
\alias{RobModel-class}
\alias{center}
\alias{center,RobModel-method}
\alias{center<-}
\alias{center<-,RobModel-method}
\alias{name,RobModel-method}
\alias{neighbor}
\alias{neighbor,RobModel-method}
\alias{neighbor<-}
\alias{neighbor<-,RobModel-method}

\title{Robust model}
\description{Class of robust models. A robust model consists 
  of  family of probability measures \code{center} and a 
  neighborhood \code{neighbor} about this family.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{center}:}{Object of class \code{"ProbFamily"} }
    \item{\code{neighbor}:}{Object of class \code{"Neighborhood"}}
  }
}
\section{Methods}{
  \describe{
    \item{center}{\code{signature(object = "RobModel")}: 
      accessor function for slot \code{center}. }
    \item{center<-}{\code{signature(object = "RobModel")}: 
      replacement function for slot \code{center}. }
    \item{neighbor}{\code{signature(object = "RobModel")}: 
      accessor function for slot \code{neighbor}. }
    \item{neighbor<-}{\code{signature(object = "RobModel")}: 
      replacement function for slot \code{neighbor}. }
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{ProbFamily-class}}, \code{\link{Neighborhood-class}}}
%\examples{}
\concept{robust model}
\keyword{classes}
\keyword{models}
\keyword{robust}