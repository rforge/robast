\name{CondTotalVarIC-class}
\docType{class}
\alias{CondTotalVarIC-class}
\alias{CallL2Fam<-,CondTotalVarIC-method}
\alias{clipLo,CondTotalVarIC-method}
\alias{clipLo<-,CondTotalVarIC-method}
\alias{clipUp,CondTotalVarIC-method}
\alias{clipUp<-,CondTotalVarIC-method}
\alias{lowerCase,CondTotalVarIC-method}
\alias{lowerCase<-,CondTotalVarIC-method}
\alias{neighborRadius,CondTotalVarIC-method}
\alias{neighborRadius<-,CondTotalVarIC-method}
\alias{neighborRadiusCurve,CondTotalVarIC-method}
\alias{neighborRadiusCurve<-,CondTotalVarIC-method}
\alias{stand,CondTotalVarIC-method}
\alias{stand<-,CondTotalVarIC-method}
\alias{generateIC,CondTotalVarNeighborhood,L2RegTypeFamily-method}
\alias{show,CondTotalVarIC-method}

\title{Conditionally centered influence curve of total variaton type}
\description{
  Class of conditionally centered (partial) influence 
  curves of contamination type for average conditional total variation
  \deqn{\eta = max(c(x), min(A x\Lambda_f, b(x)))}{
  eta = max(c(x), min(A x Lambda_f, b(x)))}
  with lower clipping function \eqn{c}, upper clipping function \eqn{b} and
  standardizing matrix \eqn{A}. \eqn{\Lambda_f}{Lambda_f} stands for
  the L2 derivative of the corresponding error distribution.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ContTotalVarIC", ...)}.
  More frequently they are created via the generating function 
  \code{ContTotalVarIC}, respectively via the method \code{generateIC}. 
}
\section{Slots}{
  \describe{
    \item{\code{CallL2Fam}:}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      regression type family. }

    \item{\code{name}:}{ object of class \code{"character"} }

    \item{\code{Curve}:}{ object of class \code{"EuclRandVarList"}}

    \item{\code{Risks}:}{ object of class \code{"list"}: 
      list of risks; cf. \code{RiskType-class}. }

    \item{\code{Infos}:}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }

    \item{\code{clipUp}:}{ object of class \code{"RealRandVariable"}:
      upper clipping function. }

    \item{\code{clipLo}:}{ object of class \code{"RealRandVariable"}:
      lower clipping function. }

    \item{\code{stand}:}{ object of class \code{"matrix"}:
      standardizing matrix. }

    \item{\code{lowerCase}:}{ object of class \code{"OptionalNumeric"}:
      optional constant for lower case solution. }

    \item{\code{neighborRadius}:}{ object of class \code{"numeric"}:
      radius of the corresponding conditional contamination neighborhood. }

    \item{\code{neighborRadiusCurve}:}{ object of class \code{"numeric"}:
      radius curve of the corresponding conditional contamination neighborhood. }
  }
}
\section{Extends}{
Class \code{"CondIC"}, directly.
Class \code{"IC"}, by class \code{"CondIC"}.
Class \code{"InfluenceCurve"}, by class \code{"CondIC"}.
}
\section{Methods}{
  \describe{
    \item{CallL2Fam<-}{\code{signature(object = "CondTotalVarIC")}: 
      replacement function for slot \code{CallL2Fam}. }

    \item{clipLo}{\code{signature(object = "CondTotalVarIC")}: 
      accessor function for slot \code{clipLo}. }

    \item{clipLo<-}{\code{signature(object = "CondTotalVarIC")}: 
      replacement function for slot \code{clipLo}. }

    \item{clipUp}{\code{signature(object = "CondTotalVarIC")}: 
      accessor function for slot \code{clipUp}. }

    \item{clipUp<-}{\code{signature(object = "CondTotalVarIC")}: 
      replacement function for slot \code{clipUp}. }

    \item{stand}{\code{signature(object = "CondTotalVarIC")}: 
      accessor function for slot \code{stand}. }

    \item{stand<-}{\code{signature(object = "CondTotalVarIC")}:
      replacement function for slot \code{stand}. }

    \item{lowerCase}{\code{signature(object = "CondTotalVarIC")}: 
      accessor function for slot \code{lowerCase}. }

    \item{lowerCase<-}{\code{signature(object = "CondTotalVarIC")}: 
      replacement function for slot \code{lowerCase}. }

    \item{neighborRadius}{\code{signature(object = "CondTotalVarIC")}: 
      accessor function for slot \code{neighborRadius}. }

    \item{neighborRadius<-}{\code{signature(object = "CondTotalVarIC")}: 
      replacement function for slot \code{neighborRadius}. }

    \item{neighborRadiusCurve}{\code{signature(object = "CondTotalVarIC")}: 
      accessor function for slot \code{neighborRadiusCurve}. }

    \item{neighborRadiusCurve<-}{\code{signature(object = "CondTotalVarIC")}: 
      replacement function for slot \code{neighborRadiusCurve}. }

    \item{generateIC}{\code{signature(neighbor = "CondTotalVarNeighborhood", L2Fam = "L2RegTypeFamily")}:
      generate an object of class \code{"CondTotalVarIC"}. Rarely called directly. }

    \item{show}{\code{signature(object = "CondTotalVarIC")} }
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}, \code{\link{CondTotalVarIC}}}
\examples{
IC1 <- new("CondTotalVarIC")
IC1
}
\concept{influence curve}
\keyword{classes}
\keyword{robust}
