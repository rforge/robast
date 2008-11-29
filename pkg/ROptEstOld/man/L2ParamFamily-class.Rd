\name{L2ParamFamily-class}
\docType{class}
\alias{L2ParamFamily-class}
\alias{FisherInfo}
\alias{FisherInfo,L2ParamFamily-method}
\alias{L2deriv}
\alias{L2deriv,L2ParamFamily-method}
\alias{L2derivSymm}
\alias{L2derivSymm,L2ParamFamily-method}
\alias{L2derivDistr}
\alias{L2derivDistr,L2ParamFamily-method}
\alias{L2derivDistrSymm}
\alias{L2derivDistrSymm,L2ParamFamily-method}
\alias{checkL2deriv,L2ParamFamily-method}
\alias{E,L2ParamFamily,EuclRandVariable,missing-method}
\alias{E,L2ParamFamily,EuclRandMatrix,missing-method}
\alias{E,L2ParamFamily,EuclRandVarList,missing-method}
\alias{plot,L2ParamFamily,ANY-method}

\title{L2 differentiable parametric family}
\description{Class of L2 differentiable parametric families.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("L2ParamFamily", ...)}.
  More frequently they are created via the generating function 
  \code{L2ParamFamily}. 
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{ object of class \code{"character"}: 
      name of the family. }
    \item{\code{distribution}:}{ object of class \code{"Distribution"}:
      member of the family. }
    \item{\code{distrSymm}:}{Object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
    \item{\code{param}:}{ object of class \code{"ParamFamParameter"}:
      parameter of the family. }
    \item{\code{props}:}{ object of class \code{"character"}:
      properties of the family. }
    \item{\code{L2deriv}:}{ object of class \code{"EuclRandVariable"}:
      L2 derivative of the family. }
    \item{\code{L2derivSymm}:}{ object of class \code{"FunSymmList"}:
      symmetry of the maps included in \code{L2deriv}. }
    \item{\code{L2derivDistr}:}{ object of class \code{"UnivarDistrList"}:
      list which includes the distribution of \code{L2deriv}. }
    \item{\code{L2derivDistrSymm}:}{ object of class \code{"DistrSymmList"}:
      symmetry of the distributions included in \code{L2derivDistr}. }
    \item{\code{FisherInfo}:}{ object of class \code{"PosDefSymmMatrix"}:
      Fisher information of the family. }
  }
}
\section{Extends}{
Class \code{"ParamFamily"}, directly.\cr
Class \code{"ProbFamily"}, by class \code{"ParamFamily"}.
}
\section{Methods}{
  \describe{
    \item{L2deriv}{\code{signature(object = "L2ParamFamily")}: 
      accessor function for \code{L2deriv}. }

    \item{L2derivSymm}{\code{signature(object = "L2ParamFamily")}: 
      accessor function for \code{L2derivSymm}. }

    \item{L2derivDistr}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{L2derivDistr}. }

    \item{L2derivDistrSymm}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{L2derivDistrSymm}. }

    \item{FisherInfo}{\code{signature(object = "L2ParamFamily")}:
      accessor function for \code{FisherInfo}. }

    \item{checkL2deriv}{\code{signature(object = "L2ParamFamily")}: 
      check centering of \code{L2deriv} and compute precision 
      of Fisher information. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandVariable", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandMatrix", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{E}{\code{signature(object = "L2ParamFamily", fun = "EuclRandVarList", cond = "missing")}:
      expectation of \code{fun} under the distribution of \code{object}. }

    \item{plot}{\code{signature(x = "L2ParamFamily")}: 
      plot of \code{distribution} and \code{L2deriv}. }
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{L2ParamFamily}}, \code{\link{ParamFamily-class}}}
\examples{
F1 <- new("L2ParamFamily")
plot(F1)
}
\concept{parametric family}
\keyword{classes}
\keyword{models}
