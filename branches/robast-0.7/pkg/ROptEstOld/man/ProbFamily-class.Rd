\name{ProbFamily-class}
\docType{class}
\alias{ProbFamily-class}
\alias{addProp<-}
\alias{addProp<-,ProbFamily-method}
\alias{distribution}
\alias{distribution,ProbFamily-method}
\alias{distrSymm}
\alias{distrSymm,ProbFamily-method}
\alias{name,ProbFamily-method}
\alias{name<-,ProbFamily-method}
\alias{props}
\alias{props,ProbFamily-method}
\alias{props<-}
\alias{props<-,ProbFamily-method}

\title{Family of probability measures}
\description{Class of families of probability measures.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}: 
      name of the family. }
    \item{\code{distribution}:}{Object of class \code{"Distribution"}: 
      member of the family. }
    \item{\code{distrSymm}:}{Object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
    \item{\code{props}:}{Object of class \code{"character"}: 
      properties of the family. }
  }
}
\section{Methods}{
  \describe{
    \item{name}{\code{signature(object = "ProbFamily")}: 
      accessor function for slot \code{name}. }

    \item{name<-}{\code{signature(object = "ProbFamily")}: 
      replacement function for slot \code{name}. }

    \item{distribution}{\code{signature(object = "ProbFamily")}: 
      accessor function for slot \code{distribution}. }

    \item{distrSymm}{\code{signature(object = "ProbFamily")}: 
      accessor function for slot \code{distrSymm}. }

    \item{props}{\code{signature(object = "ProbFamily")}: 
      accessor function for slot \code{props}. }

    \item{props<-}{\code{signature(object = "ProbFamily")}: 
      replacement function for slot \code{props}. }

    \item{addProp<-}{\code{signature(object = "ProbFamily")}: 
      add a property to slot \code{props}. }
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distr]{Distribution-class}}}
%\examples{}
\concept{family}
\concept{probability measure}
\keyword{classes}
\keyword{models}
