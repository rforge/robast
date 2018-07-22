\name{checkmakeIC-methods}
\docType{methods}
\alias{checkIC,IC,ParetoFamily-method}
\alias{checkIC,IC,GParetoFamily-method}
\alias{checkIC,IC,GEVFamily-method}
\alias{checkIC,IC,GEVFamilyMuUnknown-method}
\alias{makeIC,IC,ParetoFamily-method}
\alias{makeIC,IC,GParetoFamily-method}
\alias{makeIC,IC,GEVFamily-method}
\alias{makeIC,IC,GEVFamilyMuUnknown-method}

\title{Methods for Functions checkIC and makeIC in Package `RobExtremes' }

\description{\code{checkIC} checks accuracy of the centering
and Fisher consistency condition of an IC, \code{makeIC},
by centering and restandardizing warrants these conditions.}

\section{Methods}{\describe{
\item{checkIC}{\code{signature(IC="IC", L2Fam = "ParetoFamily")}:
      To enhance accuracy, the method for \code{"ParetoFamily"} uses
      integration via the quantile transform, i.e., \eqn{E[h(X)]}
      for a random variable \eqn{X\sim F}{X~F} with quantil function \eqn{q}
      is computed as \eqn{\int_0^1 h(q(s))\,ds}{integral(from=0, to=1, integrand=h(q(s)))}
}
\item{checkIC}{\code{signature(IC="IC", L2Fam = "GParetoFamily")}:
      As for \code{"ParetoFamily"}, to enhance accuracy,
      the method for \code{"GParetoFamily"} uses
      integration via the quantile transform.}
\item{checkIC}{\code{signature(IC="IC", L2Fam = "GEVFamily")}:
      As for \code{"ParetoFamily"}, to enhance accuracy,
      the method for \code{"GEVFamily"} uses
      integration via the quantile transform.}
\item{checkIC}{\code{signature(IC="IC", L2Fam = "GEVFamilyMuUnknown")}:
      As for \code{"ParetoFamily"}, to enhance accuracy,
      the method for \code{"GEVFamilyMuUnknown"} uses
      integration via the quantile transform.}
\item{makeIC}{\code{signature(IC="IC", L2Fam = "ParetoFamily")}:
      As with \code{"checkIC"}, to enhance accuracy,
      the method for \code{"makeIC"} for \code{"ParetoFamily"} uses
      integration via the quantile transform.}
\item{makeIC}{\code{signature(IC="IC", L2Fam = "GParetoFamily")}:
      As for \code{"ParetoFamily"}, to enhance accuracy,
      the method for \code{"GParetoFamily"} uses
      integration via the quantile transform.}
\item{makeIC}{\code{signature(IC="IC", L2Fam = "GEVFamily")}:
      As for \code{"ParetoFamily"}, to enhance accuracy,
      the method for \code{"GEVFamily"} uses
      integration via the quantile transform.}
\item{makeIC}{\code{signature(IC="IC", L2Fam = "GEVFamilyMuUnknown")}:
      As for \code{"ParetoFamily"}, to enhance accuracy,
      the method for \code{"GEVFamilyMuUnknown"} uses
      integration via the quantile transform.}
}}

\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
\seealso{\code{\link[RobAStBase]{checkIC}},\code{\link[RobAStBase:makeIC-methods]{makeIC}}}
\concept{influence curve}
\keyword{robust}
