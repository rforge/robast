\name{Av1CondContIC-class}
\docType{class}
\alias{Av1CondContIC-class}
\alias{CallL2Fam<-,Av1CondContIC-method}
\alias{cent,Av1CondContIC-method}
\alias{cent<-,Av1CondContIC-method}
\alias{clip,Av1CondContIC-method}
\alias{clip<-,Av1CondContIC-method}
\alias{lowerCase,Av1CondContIC-method}
\alias{lowerCase<-,Av1CondContIC-method}
\alias{neighborRadius,Av1CondContIC-method}
\alias{neighborRadius<-,Av1CondContIC-method}
\alias{stand,Av1CondContIC-method}
\alias{stand<-,Av1CondContIC-method}
\alias{generateIC,Av1CondContNeighborhood,L2RegTypeFamily-method}
\alias{show,Av1CondContIC-method}

\title{Conditionally centered influence curve of contamination type}
\description{
  Class of conditionally centered (partial) influence 
  curves of contamination type for average conditional contamination
  neighborhoods; i.e., influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = (A\Lambda - a)\min(1,b/|A\Lambda - a|)}{eta = (A Lambda - a)min(1, b/|A Lambda - a|)}
  with clipping bound \eqn{b}, centering function \eqn{a} and
  standardizing matrix \eqn{A}. \eqn{\Lambda}{Lambda} stands for
  the L2 derivative of the corresponding L2 differentiable 
  regression type family created via the call in the slot \code{CallL2Fam}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Av1CondContIC", ...)}.
  More frequently they are created via the generating function 
  \code{Av1CondContIC}, respectively via the method \code{generateIC}. 
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

    \item{\code{cent}:}{ object of class \code{"EuclRandVarList"}:
      centering function. }

    \item{\code{stand}:}{ object of class \code{"matrix"}:
      standardizing matrix. }

    \item{\code{lowerCase}:}{ object of class \code{"OptionalNumeric"}:
      optional constant for lower case solution. }

    \item{\code{neighborRadius}:}{ object of class \code{"numeric"}:
      radius of the corresponding average conditional contamination 
      neighborhood. }
  }
}
\section{Extends}{
Class \code{"CondIC"}, directly.\cr
Class \code{"IC"}, by class \code{"CondIC"}.\cr
Class \code{"InfluenceCurve"}, by class \code{"CondIC"}.
}
\section{Methods}{
  \describe{
    \item{CallL2Fam<-}{\code{signature(object = "Av1CondContIC")}: 
      replacement function for slot \code{CallL2Fam}. }

    \item{cent}{\code{signature(object = "Av1CondContIC")}: 
      accessor function for slot \code{cent}. }

    \item{cent<-}{\code{signature(object = "Av1CondContIC")}: 
      replacement function for slot \code{cent}. }

    \item{clip}{\code{signature(object = "Av1CondContIC")}: 
      accessor function for slot \code{clip}. }

    \item{clip<-}{\code{signature(object = "Av1CondContIC")}: 
      replacement function for slot \code{clip}. }

    \item{stand}{\code{signature(object = "Av1CondContIC")}: 
      accessor function for slot \code{stand}. }

    \item{stand<-}{\code{signature(object = "Av1CondContIC")}:
      replacement function for slot \code{stand}. }

    \item{lowerCase}{\code{signature(object = "Av1CondContIC")}: 
      accessor function for slot \code{lowerCase}. }

    \item{lowerCase<-}{\code{signature(object = "Av1CondContIC")}: 
      replacement function for slot \code{lowerCase}. }

    \item{neighborRadius}{\code{signature(object = "Av1CondContIC")}: 
      accessor function for slot \code{neighborRadius}. }

    \item{neighborRadius<-}{\code{signature(object = "Av1CondContIC")}: 
      replacement function for slot \code{neighborRadius}. }

    \item{generateIC}{\code{signature(neighbor = "Av1CondContNeighborhood", L2Fam = "L2RegTypeFamily")}:
      generate an object of class \code{"Av1CondContIC"}. Rarely called directly. }

    \item{show}{\code{signature(object = "Av1CondContIC")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{CondIC-class}}, \code{\link{Av1CondContIC}}}
\examples{
IC1 <- new("Av1CondContIC")
IC1
}
\concept{influence curve}
\keyword{classes}
\keyword{robust}
