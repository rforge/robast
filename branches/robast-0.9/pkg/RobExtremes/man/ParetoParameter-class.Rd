\name{ParetoParameter-class}
\docType{class}
\alias{ParetoParameter-class}
\alias{shape}
\alias{shape,ParetoParameter-method}
\alias{shape<-}
\alias{shape<-,ParetoParameter-method}
\alias{Min,ParetoParameter-method}
\alias{Min<-,ParetoParameter-method}

\title{Paramter of Pareto distributions}
\description{The class of the parameter of Pareto distributions.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParetoParameter", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{shape}}{ real number: shape parameter of 
      a Pareto distribution. }
    \item{\code{Min}}{ positive real number: Min
      parameter of a Pareto distribution. }
    \item{\code{name}}{ default name is 
      \dQuote{parameter of a Pareto distribution}. }
  }
}
\section{Extends}{
Class \code{"Parameter"}, directly.\cr
Class \code{"OptionalParameter"}, by class \code{"Parameter"}.
}
\section{Methods}{
  \describe{
    \item{shape}{\code{signature(object = "ParetoParameter")}: access method for
      slot \code{shape}. }
    \item{Min}{\code{signature(x = "ParetoParameter")}: access method for
      slot \code{Min}. }
    \item{shape<-}{\code{signature(object = "ParetoParameter")}: replace method for
      slot \code{shape}. }
    \item{Min<-}{\code{signature(x = "ParetoParameter")}: replace method for
      slot \code{Min}. }
  }
}
%\references{}
\author{Nataliya Horbenko \email{Nataliya.Horbenko@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{Pareto-class}}, \code{\link[distr]{Parameter-class}}}
\examples{
new("ParetoParameter")
}
\concept{Pareto distribution}
\keyword{distribution}
\concept{parameter}
\concept{S4 parameter class}
\keyword{models}
\concept{generating function}
