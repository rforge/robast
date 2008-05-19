\name{L2RegTypeFamily-class}
\docType{class}
\alias{L2RegTypeFamily-class}
\alias{ErrorL2deriv}
\alias{ErrorL2deriv,L2RegTypeFamily-method}
\alias{ErrorL2derivDistr}
\alias{ErrorL2derivDistr,L2RegTypeFamily-method}
\alias{ErrorL2derivSymm}
\alias{ErrorL2derivSymm,L2RegTypeFamily-method}
\alias{ErrorL2derivDistrSymm}
\alias{ErrorL2derivDistrSymm,L2RegTypeFamily-method}
\alias{FisherInfo,L2RegTypeFamily-method}
\alias{L2deriv,L2RegTypeFamily-method}
\alias{checkL2deriv,L2RegTypeFamily-method}
\alias{checkIC,IC,missing-method}
\alias{checkIC,IC,L2RegTypeFamily-method}
\alias{E,L2RegTypeFamily,EuclRandVariable,missing-method}
\alias{E,L2RegTypeFamily,EuclRandMatrix,missing-method}
\alias{E,L2RegTypeFamily,EuclRandVarList,missing-method}
\alias{show,L2RegTypeFamily-method}

\title{L2 differentiable parametric regression-type family}
\description{Class for L2 differentiable parametric regression-type families.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("L2RegTypeFamily", ...)}.
  More frequently they are created via the generating function 
  \code{L2RegTypeFamily}. 
}
\section{Slots}{
  \describe{
    \item{\code{L2deriv}:}{Object of class \code{"EuclRandVarList"}: L2 derivative. }
    \item{\code{ErrorL2deriv}:}{Object of class \code{"EuclRandVarList"}: 
      L2 derivative of \code{ErrorDistr}. }
    \item{\code{ErrorL2derivSymm}:}{Object of class \code{"FunSymmList"}:
      symmetry of \code{ErrorL2deriv}. }
    \item{\code{ErrorL2derivDistr}:}{Object of class \code{"DistrList"}: 
       distribution of \code{ErrorL2deriv}. }
    \item{\code{ErrorL2derivDistrSymm}:}{Object of class \code{"DistrSymmList"}:
      symmetry of \code{ErrorL2derivDistr}. }
    \item{\code{FisherInfo}:}{Object of class \code{"PosDefSymmMatrix"}:
      Fisher information. }
    \item{\code{ErrorDistr}:}{Object of class \code{"Distribution"}: 
      error distribution. }
    \item{\code{ErrorSymm}:}{Object of class \code{"DistributionSymmetry"}:
      symmetry of \code{ErrorDistr}. }
    \item{\code{RegDistr}:}{Object of class \code{"Distribution"}:
      regressor distribution. }
    \item{\code{RegSymm}:}{Object of class \code{"DistributionSymmetry"}:
      symmetry of \code{RegDistr}. }
    \item{\code{Regressor}:}{Object of class \code{"EuclRandVariable"}:
      regressor. }
    \item{\code{param}:}{Object of class \code{"ParamFamParameter"}:
      parameter of the family. }
    \item{\code{props}:}{Object of class \code{"character"}:
      properties of the family. }
    \item{\code{name}:}{Object of class \code{"character"}:
      name of the family. }
    \item{\code{distribution}:}{Object of class \code{"CondDistribution"}:
      conditional distribution given the regressor. }
    \item{\code{distrSymm}:}{Object of class \code{"DistributionSymmetry"}:
      symmetry of \code{distribution}. }
  }
}
\section{Extends}{
Class \code{"RegTypeFamily"}, directly.\cr
Class \code{"ParamFamily"}, by class \code{"RegTypeFamily"}.\cr
Class \code{"ProbFamily"}, by class \code{"RegTypeFamily"}.
}
\section{Methods}{
  \describe{
    \item{L2deriv}{\code{signature(object = "L2RegTypeFamily")}:
      acessor function for slot \code{L2deriv} }

    \item{FisherInfo}{\code{signature(object = "L2RegTypeFamily")}:
      acessor function for slot \code{FisherInfo} }

    \item{ErrorL2deriv}{\code{signature(object = "L2RegTypeFamily")}: 
      acessor function for slot \code{ErrorL2deriv}. }

    \item{ErrorL2derivDistr}{\code{signature(object = "L2RegTypeFamily")}: 
      acessor function for slot \code{ErrorL2derivDistr}. }

    \item{ErrorL2derivSymm}{\code{signature(object = "L2RegTypeFamily")}: 
      acessor function for slot \code{ErrorL2derivSymm} }

    \item{ErrorL2derivDistrSymm}{\code{signature(object = "L2RegTypeFamily")}: 
      acessor function for slot \code{ErrorL2derivDistrSymm} }

    \item{checkL2deriv}{\code{signature(object = "L2RegTypeFamily")}: 
      check centering of \code{L2deriv} and compute precision 
      of Fisher information. }

    \item{checkIC}{\code{signature(IC = "IC", L2Fam = "missing")}: 
      check centering and Fisher consistency of \code{IC} assuming the 
      L2-differentiable regression-type family which can be created via 
      the slot \code{CallL2Fam} of \code{IC}. }

    \item{checkIC}{\code{signature(IC = "IC", L2Fam = "L2RegTypeFamily")}: 
      check centering and Fisher consistency of \code{IC} assuming the 
      L2-differentiable regression-type family \code{L2Fam}. }

    \item{E}{\code{signature(object = "L2RegTypeFamily", fun = "EuclRandVariable", cond = "missing")}:
      expectation of \code{fun} under \code{object}. }

    \item{E}{\code{signature(object = "L2RegTypeFamily", fun = "EuclRandMatrix", cond = "missing")}:
      expectation of \code{fun} under \code{object}. }

    \item{E}{\code{signature(object = "L2RegTypeFamily", fun = "EuclRandVarList", cond = "missing")}:
      expectation of \code{fun} under \code{object}. }

    \item{show}{\code{signature(object = "L2RegTypeFamily")}}
  }
}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{RegTypeFamily-class}}}
\examples{
new("L2RegTypeFamily")
}
\concept{regression}
\keyword{classes}
\keyword{models}
