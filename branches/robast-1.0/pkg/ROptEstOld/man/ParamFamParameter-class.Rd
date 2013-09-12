\name{ParamFamParameter-class}
\docType{class}
\alias{ParamFamParameter-class}
\alias{length,ParamFamParameter-method}
\alias{main}
\alias{main,ParamFamParameter-method}
\alias{main<-}
\alias{main<-,ParamFamParameter-method}
\alias{nuisance}
\alias{nuisance,ParamFamParameter-method}
\alias{nuisance<-}
\alias{nuisance<-,ParamFamParameter-method}
\alias{show,ParamFamParameter-method}
\alias{trafo}
\alias{trafo,ParamFamParameter-method}
\alias{trafo<-}
\alias{trafo<-,ParamFamParameter-method}

\title{Parameter of a parametric family of probability measures}
\description{Class of the parameter of parametric families 
  of probability measures.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ParamFamParameter", ...)}.
  More frequently they are created via the generating function 
  \code{ParamFamParameter}. 
}
\section{Slots}{
  \describe{
    \item{\code{main}:}{Object of class \code{"numeric"}: main parameter. }
    \item{\code{nuisance}:}{Object of class \code{"OptionalNumeric"}:
      optional nuisance parameter. }
    \item{\code{trafo}:}{Object of class \code{"matrix"}: 
      transformation of the parameter. }
    \item{\code{name}:}{Object of class \code{"character"}:
      name of the parameter. }
  }
}
\section{Extends}{
Class \code{"Parameter"}, directly.\cr
Class \code{"OptionalParameter"}, by class \code{"Parameter"}.
}
\section{Methods}{
  \describe{
    \item{main}{\code{signature(object = "ParamFamParameter")}: 
      accessor function for slot \code{main}. }

    \item{main<-}{\code{signature(object = "ParamFamParameter")}: 
      replacement function for slot \code{main}. }

    \item{nuisance}{\code{signature(object = "ParamFamParameter")}: 
      accessor function for slot \code{nuisance}. }

    \item{nuisance<-}{\code{signature(object = "ParamFamParameter")}: 
      replacement function for slot \code{nuisance}. }

    \item{trafo}{\code{signature(object = "ParamFamParameter")}: 
      accessor function for slot \code{trafo}. }

    \item{trafo<-}{\code{signature(object = "ParamFamParameter")}: 
      replacement function for slot \code{trafo}. }

    \item{length}{\code{signature(x = "ParamFamParameter")}: 
      sum of the lengths of \code{main} and \code{nuisance}. }

    \item{show}{\code{signature(object = "ParamFamParameter")}}
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distr]{Parameter-class}}}
\examples{
new("ParamFamParameter")
}
\concept{parameter}
\keyword{classes}
