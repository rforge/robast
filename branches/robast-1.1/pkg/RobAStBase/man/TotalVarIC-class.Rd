\name{TotalVarIC-class}
\docType{class}
\alias{TotalVarIC-class}
\alias{CallL2Fam<-,TotalVarIC-method}
\alias{clipLo}
\alias{clip,TotalVarIC-method}
\alias{clipLo,TotalVarIC-method}
\alias{clipLo<-}
\alias{clipLo<-,TotalVarIC-method}
\alias{clipUp}
\alias{clipUp,TotalVarIC-method}
\alias{clipUp<-}
\alias{clipUp<-,TotalVarIC-method}
\alias{lowerCase<-,TotalVarIC-method}
\alias{neighbor,TotalVarIC-method}
\alias{show,TotalVarIC-method}
\alias{stand<-,TotalVarIC-method}
\alias{generateIC,TotalVarNeighborhood,L2ParamFamily-method}

\title{Influence curve of total variation type}
\description{Class of (partial) influence curves of total variation type.
  i.e., an influence curves \eqn{\eta}{eta} of the form
  \deqn{\eta = c \vee A\Lambda \wedge d}{eta = max(c, min(A Lambda, d))}
  with lower clipping bound \eqn{c}, upper clipping bound \eqn{d} and
  standardizing matrix \eqn{A}. \eqn{\Lambda}{Lambda} stands for
  the L2 derivative of the corresponding L2 differentiable 
  parametric family which can be created via \code{CallL2Fam}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("TotalVarIC", ...)}.
  More frequently they are created via the generating function 
  \code{TotalVarIC}, respectively via the method \code{generateIC}. 
}
\section{Slots}{
  \describe{
    \item{\code{CallL2Fam}}{ object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      parametric family. }

    \item{\code{name}}{ object of class \code{"character"}. }

    \item{\code{Curve}}{ object of class \code{"EuclRandVarList"}.}

    \item{\code{modifyIC}}{ Object of class \code{"OptionalFunction"}:
      function of two arguments, which are an L2 parametric family 
      and an optional influence curve. Returns an object of 
      class \code{"IC"}. This slot is mainly used for internal
      computations! }

    \item{\code{Risks}}{ object of class \code{"list"}: 
      list of risks; cf. \code{\link[distrMod]{RiskType-class}}. }

    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }

    \item{\code{clipLo}}{ object of class \code{"numeric"}: 
      lower clipping bound. }

    \item{\code{clipUp}}{ object of class \code{"numeric"}: 
      upper clipping bound. }

    \item{\code{stand}}{ object of class \code{"matrix"}:
      standardizing matrix. }

    \item{\code{weight}}{ object of class \code{"BdStWeight"}:
      weight function }

    \item{\code{biastype}}{ object of class \code{"BiasType"}:
      bias type (symmetric/onsided/asymmetric) }

    \item{\code{normtype}}{ object of class \code{"NormType"}:
      norm type (Euclidean, information/self-standardized)}

    \item{\code{neighborRadius}}{ object of class \code{"numeric"}:
      radius of the corresponding (unconditional) contamination 
      neighborhood. }
  }
}
\section{Extends}{
Class \code{"HampIC"}, directly.\cr
Class \code{"IC"}, by class \code{"HampIC"}.\cr
Class \code{"InfluenceCurve"}, by class \code{"IC"}.
}
\section{Methods}{
  \describe{
    \item{CallL2Fam<-}{\code{signature(object = "TotalVarIC")}: 
      replacement function for slot \code{CallL2Fam}. }

    \item{clipLo}{\code{signature(object = "TotalVarIC")}: 
      accessor function for slot \code{clipLo}. }

    \item{clipLo<-}{\code{signature(object = "TotalVarIC")}: 
      replacement function for slot \code{clipLo}. }

    \item{clipUp}{\code{signature(object = "TotalVarIC")}: 
      accessor function for slot \code{clipUp}. }

    \item{clipUp<-}{\code{signature(object = "TotalVarIC")}: 
      replacement function for slot \code{clipUp}. }

    \item{clip}{\code{signature(x1 = "TotalVarIC")}:
      returns \code{clipUp-clipLo}. }

    \item{stand<-}{\code{signature(object = "TotalVarIC")}:
      replacement function for slot \code{stand}. }

    \item{lowerCase<-}{\code{signature(object = "TotalVarIC")}:
      replacement function for slot \code{lowerCase}. }

    \item{neighbor}{\code{signature(object = "TotalVarIC")}: 
      generates an object of class \code{"TotalVarNeighborhood"} with
      radius given in slot \code{neighborRadius}. }

    \item{generateIC}{\code{signature(neighbor = "TotalVarNeighborhood", L2Fam = "L2ParamFamily")}:
      generate an object of class \code{"TotalVarIC"}. Rarely called directly. }

    \item{show}{\code{signature(object = "TotalVarIC")}}
  }
}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{IC-class}}, \code{\link{ContIC}}, \code{\link{HampIC-class}}}
\examples{
IC1 <- new("TotalVarIC")
plot(IC1)
}
\concept{influence curve}
\keyword{classes}
\keyword{robust}
