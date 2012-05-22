\name{interpolRisk-class}
\docType{class}
\alias{interpolRisk-class}
\alias{OMSRRisk-class}
\alias{RMXRRisk-class}
\alias{MBRRisk-class}
\alias{OMSRRisk}
\alias{RMXRRisk}
\alias{MBRRisk}

\title{Interpolated Risks}
\description{Class of risks for which algorithms dispatch to speed-up algorithms}
\section{Objects from the Class}{\code{interpolRisk} is a virtual class: No objects may be created from it.
the other classes are generated via generating functions.
}
\section{Slots}{
  \describe{
    \item{\code{type}}{Object of class \code{"character"}:
      type of risk. (Inherited from \code{RiskType}).}
  }
}
\details{
The main purpose of classes \code{OMSRRisk}, \code{MBRRisk}, and
\code{RMXRRisk} is to help to dispatch into speed-up algorithms later
in function \code{roptest}. In all these risks, we assume
convex contamination neighborhoods.
\code{OMSRRisk} stands for optimal MSE-robust estimation (where we assume
a radius r of \code{0.5}), \code{RMXRRisk} stands for optimal
optimally RMX-robust estimation and \code{MBRRisk} stands for optimal
Bias-robust estimation.
}
%\references{}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
\examples{
new("OMSRRisk")
OMSRRisk()
RMXRRisk()
MBRRisk()
}

%\note{}
%\seealso{}
%\examples{}
\concept{risk}
\keyword{classes}
