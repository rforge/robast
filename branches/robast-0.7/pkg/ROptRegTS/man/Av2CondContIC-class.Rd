\name{Av2CondContIC-class}
\docType{class}
\alias{Av2CondContIC-class}
\alias{CallL2Fam<-,Av2CondContIC-method}
\alias{cent,Av2CondContIC-method}
\alias{cent<-,Av2CondContIC-method}
\alias{clip,Av2CondContIC-method}
\alias{clip<-,Av2CondContIC-method}
\alias{lowerCase,Av2CondContIC-method}
\alias{lowerCase<-,Av2CondContIC-method}
\alias{neighborRadius,Av2CondContIC-method}
\alias{neighborRadius<-,Av2CondContIC-method}
\alias{stand,Av2CondContIC-method}
\alias{stand<-,Av2CondContIC-method}
\alias{generateIC,Av2CondContNeighborhood,L2RegTypeFamily-method}
\alias{show,Av2CondContIC-method}

\title{Conditionally centered influence curve of contamination type}
\description{Class of conditionally centered (partial) influence 
  curves of contamination type for average square conditional contamination
  neighborhoods; i.e., influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = A K^{-1} x(\Lambda_f - z)\min(1,c/|\Lambda_f - z|)}{eta = A K^(-1) x(Lambda_f - z)min(1, c/|Lambda_f - a|)}
  with \eqn{K=E xx^\tau}{K=E x t(x)}, clipping bound \eqn{c}, centering constant 
  \eqn{z} and standardizing constant \eqn{A}. \eqn{\Lambda_f}{Lambda_f} stands for
  the L2 derivative of the corresponding error distribution.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Av2CondContIC", ...)}.
  More frequently they are created via the generating function 
  \code{Av2CondContIC}, respectively via the method \code{generateIC}.
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

    \item{\code{clip}:}{ object of class \code{"numeric"}: 
      clipping bound. }

    \item{\code{cent}:}{ object of class \code{"numeric"}:
      centering constant. }

    \item{\code{stand}:}{ object of class \code{"numeric"}:
      standardizing constant. }

    \item{\code{lowerCase}:}{ object of class \code{"OptionalNumeric"}:
      optional constant for lower case solution. }

    \item{\code{neighborRadius}:}{ object of class \code{"numeric"}:
      radius of the corresponding average conditional contamination 
      neighborhood. }
  }
}
\section{Extends}{
Class \code{"CondIC"}, directly.
Class \code{"IC"}, by class \code{"CondIC"}.
Class \code{"InfluenceCurve"}, by class \code{"CondIC"}.
}
\section{Methods}{
  \describe{
    \item{CallL2Fam<-}{\code{signature(object = "Av2CondContIC")}: 
      replacement function for slot \code{CallL2Fam}. }

    \item{cent}{\code{signature(object = "Av2CondContIC")}: 
      accessor function for slot \code{cent}. }

    \item{cent<-}{\code{signature(object = "Av2CondContIC")}: 
      replacement function for slot \code{cent}. }

    \item{clip}{\code{signature(object = "Av2CondContIC")}: 
      accessor function for slot \code{clip}. }

    \item{clip<-}{\code{signature(object = "Av2CondContIC")}: 
      replacement function for slot \code{clip}. }

    \item{stand}{\code{signature(object = "Av2CondContIC")}: 
      accessor function for slot \code{stand}. }

    \item{stand<-}{\code{signature(object = "Av2CondContIC")}:
      replacement function for slot \code{stand}. }

    \item{lowerCase}{\code{signature(object = "Av2CondContIC")}: 
      accessor function for slot \code{lowerCase}. }

    \item{lowerCase<-}{\code{signature(object = "Av2CondContIC")}: 
      replacement function for slot \code{lowerCase}. }

    \item{neighborRadius}{\code{signature(object = "Av2CondContIC")}: 
      accessor function for slot \code{neighborRadius}. }

    \item{neighborRadius<-}{\code{signature(object = "Av2CondContIC")}: 
      replacement function for slot \code{neighborRadius}. }

    \item{generateIC}{\code{signature(neighbor = "Av2CondContNeighborhood", L2Fam = "L2RegTypeFamily")}:
      generate an object of class \code{"Av2CondContIC"}. Rarely called directly. }

    \item{show}{\code{signature(object = "Av2CondContIC")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}, \code{\link{Av2CondContIC}}}
\examples{
IC1 <- new("Av2CondContIC")
IC1
}
\concept{influence curve}
\keyword{classes}
\keyword{robust}
