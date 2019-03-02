\name{getRiskFctBV-methods}
\docType{methods}
\alias{getRiskFctBV}
\alias{getRiskFctBV-methods}
\alias{getRiskFctBV,asGRisk,ANY-method}
\alias{getRiskFctBV,asMSE,ANY-method}
\alias{getRiskFctBV,interpolRisk,ANY-method}
\alias{getRiskFctBV,asSemivar,onesidedBias-method}
\alias{getRiskFctBV,asSemivar,asymmetricBias-method}
\title{Methods for Function getRiskFctBV  in Package `RobAStBase'}
\description{getRiskFctBV for a given object of S4 class \code{asGRisk}
       returns a function in bias and variance to compute the asymptotic
       risk.}

\section{Methods}{\describe{

\item{getRiskFctBV}{\code{signature(risk = "asGRisk", biastype = "ANY")}:
  returns an error that the respective method is not yet implemented. }
\item{getRiskFctBV}{\code{signature(risk = "asMSE", biastype = "ANY")}:
  returns a function with arguments \code{bias} and \code{variance}
  to compute the asymptotic MSE for a given ALE at a situation where it has bias \code{bias}
  (including the radius!) and variance \code{variance}. }
\item{getRiskFctBV}{\code{signature(risk = "asSemivar", biastype = "onesidedBias")}:
  returns a function with arguments \code{bias} and \code{variance}
  to compute the asymptotic semivariance error, i.e. \eqn{E[(S_n-\theta)_+^2]}{E[max(Sn-theta,0)^2]}
  resp. \eqn{E[(S_n-\theta)_-^2]}{E[max(theta-Sn,0)^2]}, for a given ALE \eqn{S_n}{Sn}
  at a  situation where it has one-sided bias \code{bias} (including the radius!)
  and variance \code{variance}. }
\item{getRiskFctBV}{\code{signature(risk = "asSemivar", biastype = "asymmetricBias")}:
  returns a function with arguments \code{bias} and \code{variance}
  to compute the asymptotic semivariance error, i.e.
  \eqn{E[\nu_1 (S_n-\theta)_+^2+\nu_2(S_n-\theta)_-^2]}{nu1*max(Sn-theta,0)^2+nu2*max(theta-Sn,0)^2]}
  for a given ALE \eqn{S_n}{Sn} at a situation where it has one-sided bias
  \code{bias} (including the radius!) and variance \code{variance}. }
}}
\examples{
myrisk <- asMSE()
getRiskFctBV(myrisk)
}
\concept{risk}
\keyword{classes}


