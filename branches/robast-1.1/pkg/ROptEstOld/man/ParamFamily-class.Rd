\name{ParamFamily-class}
\docType{class}
\alias{ParamFamily-class}
\alias{main,ParamFamily-method}
\alias{nuisance,ParamFamily-method}
\alias{param,ParamFamily-method}
\alias{plot,ParamFamily,ANY-method}
\alias{show,ParamFamily-method}
\alias{trafo,ParamFamily-method}

\title{Parametric family of probability measures.}
\description{Class of parametric families of probability measures.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ParamFamily", ...)}.
  More frequently they are created via the generating function 
  \code{ParamFamily}. 
}
\section{Slots}{
  \describe{
    \item{\code{param}:}{Object of class \code{"ParamFamParameter"}:
      parameter of the family. }
    \item{\code{name}:}{Object of class \code{"character"}:
      name of the family. }
    \item{\code{distribution}:}{Object of class \code{"Distribution"}:
      member of the family.}
    \item{\code{distrSymm}:}{Object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
    \item{\code{props}:}{Object of class \code{"character"}:
      properties of the family. }
  }
}
\section{Extends}{
Class \code{"ProbFamily"}, directly.
}
\section{Methods}{
  \describe{
    \item{main}{\code{signature(object = "ParamFamily")}: 
      wrapped accessor function for slot \code{main} of
      slot \code{param}. }

    \item{nuisance}{\code{signature(object = "ParamFamily")}: 
      wrapped accessor function for slot \code{nuisance} 
      of slot \code{param}. }

    \item{trafo}{\code{signature(object = "ParamFamily")}: 
      wrapped accessor function for slot \code{trafo} 
      of slot \code{param}. }

    \item{param}{\code{signature(object = "ParamFamily")}: 
      accessor function for slot \code{param}. }


    \item{plot}{\code{signature(x = "ParamFamily")}: 
      plot of slot \code{distribution}. }

    \item{show}{\code{signature(object = "ParamFamily")}}
  }
}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distr]{Distribution-class}}}
\examples{
F1 <- new("ParamFamily") # prototype
plot(F1)
}
\concept{parameteric family}
\keyword{classes}
\keyword{models}
