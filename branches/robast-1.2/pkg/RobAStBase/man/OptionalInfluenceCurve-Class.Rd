\name{OptionalInfluenceCurve-class}
\docType{class}
\alias{OptionalInfluenceCurve-class}
\alias{OptionalInfluenceCurveOrCall-class}
\alias{OptionalpICList-class}
\alias{StartClass-class}
\alias{pICList-class}
\alias{show,pICList-method}
\alias{show,OptionalpICList-method}

\title{Some helper Classes in package 'RobAStBase'}
\description{Some helper Classes in package 'RobAStBase':
Classes \code{OptionalInfluenceCurve}, \code{OptionalpICList},
\code{StartClass}, \code{pICList}}
\section{Class Unions}{
  \code{OptionalInfluenceCurve} is a class union of classes
  \code{InfluenceCurve} and \code{NULL};
  \code{OptionalInfluenceCurveOrCall} is a class union of classes
  \code{InfluenceCurve}, \code{call}, and \code{NULL} --- it is the slot
  class of slot \code{pIC} in \code{ALEstimate};
  \code{OptionalpICList} is a class union of classes
  \code{pICList} and \code{NULL} --- it is the slot
  class of slot \code{pICList} in \code{kStepEstimate};
  \code{StartClass} is a class union of classes
  \code{function}, \code{numeric} and \code{Estimate} --- it is the slot
  class of slot \code{start} in \code{kStepEstimate}.
}
\section{List Classes}{
 \code{pICList} is a descendant of class \code{list} which requires
 its members ---if any--- to be of class \code{pIC}.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "OptionalpICList")}:
      particular show-method. }
    \item{show}{\code{signature(object = "pICList")}:
      particular show-method. }
}}
\references{
  Hampel et al. (1986) \emph{Robust Statistics}.
  The Approach Based on Influence Functions. New York: Wiley.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{InfluenceCurve}}, \code{\link[distrMod]{RiskType-class}}}
\concept{influence curve}
\keyword{classes}
\keyword{robust}
