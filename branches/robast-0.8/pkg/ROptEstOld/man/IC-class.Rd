\name{IC-class}
\docType{class}
\alias{IC-class}
\alias{CallL2Fam}
\alias{CallL2Fam,IC-method}
\alias{CallL2Fam<-}
\alias{CallL2Fam<-,IC-method}
\alias{checkIC,IC,missing-method}
\alias{checkIC,IC,L2ParamFamily-method}
\alias{evalIC,IC,numeric-method}
\alias{evalIC,IC,matrix-method}
\alias{infoPlot,IC-method}
\alias{plot,IC,ANY-method}
\alias{show,IC-method}

\title{Influence curve}
\description{Class of (partial) influence curves.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("IC", ...)}.
  More frequently they are created via the generating function 
  \code{IC}. 
}
\section{Slots}{
  \describe{
    \item{\code{CallL2Fam}:}{Object of class \code{"call"}:
      creates an object of the underlying L2-differentiable 
      parametric family. }
    \item{\code{name}:}{Object of class \code{"character"}. }
    \item{\code{Curve}:}{Object of class \code{"EuclRandVarList"}.}
    \item{\code{Risks}:}{Object of class \code{"list"}: 
      list of risks; cf. \code{\link{RiskType-class}}. }
    \item{\code{Infos}:}{Object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
  }
}
\section{Extends}{
Class \code{"InfluenceCurve"}, directly.
}
\section{Methods}{
  \describe{
    \item{CallL2Fam}{\code{signature(object = "IC")}: 
      accessor function for slot \code{CallL2Fam}. }

    \item{CallL2Fam<-}{\code{signature(object = "IC")}: 
      replacement function for slot \code{CallL2Fam}. }

    \item{checkIC}{\code{signature(IC = "IC", L2Fam = "missing")}: 
      check centering and Fisher consistency of \code{IC} assuming
      the L2-differentiable parametric family which can 
      be generated via the slot \code{CallL2Fam} of \code{IC}. }

    \item{checkIC}{\code{signature(IC = "IC", L2Fam = "L2ParamFamily")}: 
      check centering and Fisher consistency of \code{IC} assuming
      the L2-differentiable parametric family \code{L2Fam}. }

    \item{evalIC}{\code{signature(IC = "IC", x = "numeric")}: 
      evaluate \code{IC} at \code{x}. }

    \item{evalIC}{\code{signature(IC = "IC", x = "matrix")}: 
      evaluate \code{IC} at the rows of \code{x}. }

    \item{infoPlot}{\code{signature(object = "IC")}: 
      Plot absolute and relative information of \code{IC}. }

    \item{plot}{\code{signature(x = "IC")}}

    \item{show}{\code{signature(object = "IC")}}
  }
}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}. 
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve-class}}, \code{\link{IC}}}
\examples{
IC1 <- new("IC")
plot(IC1)
}
\concept{influence curve}
\keyword{classes}
