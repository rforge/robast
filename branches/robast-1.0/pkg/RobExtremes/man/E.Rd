\name{E}
\alias{E}
\alias{DistributionsIntegratingByQuantiles-class}
\alias{E-methods}
\alias{E,DistributionsIntegratingByQuantiles,function,missing-method}
\alias{E,Gumbel,missing,missing-method}
\alias{E,GPareto,missing,missing-method}
\alias{E,GPareto,function,missing-method}
\alias{E,GEV,function,missing-method}
\alias{E,Weibull,function,missing-method}
\alias{E,GEV,missing,missing-method}
\alias{E,Pareto,missing,missing-method}
\alias{E,Gammad,function,missing-method}
\alias{E,Pareto,function,missing-method}

\title{Generic Function for the Computation of (Conditional) Expectations}
\description{
  Generic function for the computation of (conditional) expectations.
}
\usage{
E(object, fun, cond, ...)

\S4method{E}{GEV,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{DistributionsIntegratingByQuantiles,function,missing}(object,
         fun, low = NULL, upp = NULL,
         rel.tol= getdistrExOption("ErelativeTolerance"),
         lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
         upperTruncQuantile = getdistrExOption("EupperTruncQuantile"),
         IQR.fac = max(1e4,getdistrExOption("IQR.fac")), ...)
\S4method{E}{Gumbel,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{GPareto,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{GPareto,function,missing}(object, fun, low = NULL, upp = NULL,
             rel.tol= getdistrExOption("ErelativeTolerance"),
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"),
             IQR.fac = max(1e4,getdistrExOption("IQR.fac")), ...)
\S4method{E}{Pareto,missing,missing}(object, low = NULL, upp = NULL, ...)
}
\arguments{
  \item{object}{ object of class \code{"Distribution"}}
  \item{fun}{ if missing the (conditional) expectation is computed
    else the (conditional) expection of \code{fun} is computed. }
  \item{cond}{ if not missing the conditional expectation 
    given \code{cond} is computed. }
  \item{rel.tol}{relative tolerance for \code{distrExIntegrate}.}
  \item{low}{lower bound of integration range.}
  \item{upp}{upper bound of integration range.}
  \item{lowerTruncQuantile}{lower quantile for quantile based integration range.}
  \item{upperTruncQuantile}{upper quantile for quantile based integration range.}
  \item{IQR.fac}{factor for scale based integration range (i.e.; 
  median of the distribution \eqn{\pm}{+-}\code{IQR.fac}\eqn{\times}{*}IQR).}
  \item{\dots}{ additional arguments to \code{fun} }
}
\details{The precision of the computations can be controlled via 
  certain global options; cf. \code{\link{distrExOptions}}. 
  Also note that arguments \code{low} and \code{upp} should be given as
  named arguments in order to prevent them to be matched by arguments
  \code{fun} or \code{cond}. Also the result, when arguments 
  \code{low} or \code{upp} is given, is the \emph{unconditional value} of the
  expectation; no conditioning with respect to \code{low <= object <= upp}
  is done. To be able to use integration after transformation via the
  respective probability transformation to [0,1], we introduce a class union
  \code{"DistributionsIntegratingByQuantiles"}, which currently comprises
  classes \code{"GPareto"}, \code{"Pareto"}, \code{"Weibull"}, \code{"GEV"}.
  In addition, the specific method for \code{"GPareto", "function", "missing"}
  uses integration on [0,1] via the substitution method (y := log(x)).
  }

\value{
  The expectation is computed.
}
\section{Methods}{
\describe{
  \item{object = "Gumbel", fun = "missing", cond = "missing":}{
    exact evaluation using explicit expressions.}
  \item{object = "GPareto", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "DistributionsIntegratingByQuantiles",
        fun = "function", cond = "missing":}{
    use probability transform, i.e., a substitution \code{y = p(object)(x)}
    for numerical integration.}
  \item{object = "GPareto", fun = "function", cond = "missing":}{
    use substitution method (y := log(x)) for numerical integration.}
  \item{object = "Pareto", fun = "missing", cond = "missing":}{
    exact evaluation using explicit expressions.}
}}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de} and Peter Ruckdeschel \email{peter.ruckdeschel@uni-bayreuth.de}}
\seealso{\code{\link{distrExIntegrate}}, \code{\link{m1df}}, \code{\link{m2df}},
    \code{\link[distr]{Distribution-class}}}
\examples{
GP <- GPareto(shape=0.3)

E(GP)
E(GP, fun = function(x){2*x^2}) ## uses the log trafo

P <- Pareto()
E(P)
E(P,fun = function(x){1/(x^2+1)})

}
\concept{expectation}
\concept{integration}
\keyword{methods}
\keyword{distribution}
\concept{functional}
\concept{moment}
