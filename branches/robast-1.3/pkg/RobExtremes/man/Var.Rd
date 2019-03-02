\name{var}
\alias{var}
\alias{var-methods}
\alias{var,Gumbel-method}
\alias{var,GPareto-method}
\alias{var,GEV-method}
\alias{var,Pareto-method}
\alias{median}
\alias{median-methods}
\alias{median,Gumbel-method}
\alias{median,GEV-method}
\alias{median,GPareto-method}
\alias{median,Pareto-method}
\alias{IQR}
\alias{IQR-methods}
\alias{IQR,Gumbel-method}
\alias{IQR,GPareto-method}
\alias{IQR,GEV-method}
\alias{IQR,Pareto-method}
\alias{skewness}
\alias{skewness-methods}
\alias{skewness,Gumbel-method}
\alias{skewness,GEV-method}
\alias{skewness,GPareto-method}
\alias{skewness,Pareto-method}
\alias{kurtosis}
\alias{kurtosis-methods}
\alias{kurtosis,Gumbel-method}
\alias{kurtosis,GEV-method}
\alias{kurtosis,GPareto-method}
\alias{kurtosis,Pareto-method}
\alias{Sn}
\alias{Sn-methods}
\alias{Sn,ANY-method}
\alias{Sn,UnivariateDistribution-method}
\alias{Sn,Norm-method}
\alias{Sn,AffLinDistribution-method}
\alias{Sn,GPareto-method}
\alias{Sn,Gammad-method}
\alias{Sn,Weibull-method}
\alias{Sn,GEV-method}
\alias{Sn,Pareto-method}
\alias{Sn,DiscreteDistribution-method}
\alias{Qn}
\alias{Qn-methods}
\alias{Qn,ANY-method}
\alias{Qn,UnivariateDistribution-method}
\alias{Qn,Norm-method}
\alias{Qn,DiscreteDistribution-method}
\alias{Qn,AffLinDistribution-method}

\title{Generic Functions for the Computation of Functionals}
\description{
  Generic functions for the computation of functionals on distributions.
}
\usage{
IQR(x, ...)

\S4method{IQR}{Gumbel}(x)
\S4method{IQR}{GEV}(x)
\S4method{IQR}{GPareto}(x)
\S4method{IQR}{Pareto}(x)

median(x, ...)

\S4method{median}{Gumbel}(x)
\S4method{median}{GEV}(x)
\S4method{median}{GPareto}(x)
\S4method{median}{Pareto}(x)

var(x, ...)

\S4method{var}{Gumbel}(x, ...)
\S4method{var}{GEV}(x, ...)
\S4method{var}{GPareto}(x, ...)
\S4method{var}{Pareto}(x, ...)

skewness(x, ...)
\S4method{skewness}{Gumbel}(x, ...)
\S4method{skewness}{GEV}(x, ...)
\S4method{skewness}{GPareto}(x, ...)
\S4method{skewness}{Pareto}(x, ...)

kurtosis(x, ...)
\S4method{kurtosis}{Gumbel}(x, ...)
\S4method{kurtosis}{GEV}(x, ...)
\S4method{kurtosis}{GPareto}(x, ...)
\S4method{kurtosis}{Pareto}(x, ...)

Sn(x, ...)
\S4method{Sn}{ANY}(x,  ...)
\S4method{Sn}{UnivariateDistribution}(x, low = 0, upp = 1.01, accuracy = 1000, ...)
\S4method{Sn}{DiscreteDistribution}(x,  ...)
\S4method{Sn}{AffLinDistribution}(x,  ...)
\S4method{Sn}{Norm}(x,  ...)
\S4method{Sn}{GPareto}(x,  ...)
\S4method{Sn}{Pareto}(x,  ...)
\S4method{Sn}{GEV}(x,  ...)
\S4method{Sn}{Gammad}(x,  ...)
\S4method{Sn}{Weibull}(x,  ...)

Qn(x, ...)
\S4method{Qn}{ANY}(x,  ...)
\S4method{Qn}{UnivariateDistribution}(x, q00 = NULL, ...)
\S4method{Qn}{AffLinDistribution}(x, ...)
\S4method{Qn}{DiscreteDistribution}(x,  ...)
\S4method{Qn}{Norm}(x,  ...)
}

\arguments{
  \item{x}{ object of class \code{"UnivariateDistribution"}}
  \item{\dots}{ additional arguments to \code{fun} or \code{E}}
  \item{q00}{numeric or NULL: determines search interval (from \code{-q00} to \code{q00})
             for \code{Qn}; if \code{NULL} (default) \code{q00}
             is set to \code{10*q(x)(3/4)} internally. }
  \item{low}{numeric; lower bound for search interval for median(abs(x-Y)) where
             Y (a real constant) runs over the range of x; defaults to \code{0}. }
  \item{upp}{numeric; upper bound for search interval for median(abs(x-Y)) where
             Y (a real constant) runs over the range of x; defaults to \code{1.01}.
            Is used internally as \code{upp*(mad(x)+abs(median(x)-Y))}. }
  \item{accuracy}{numeric; number of grid points for \code{Sn}; defaults to \code{1000}. }
  }
\value{
  The value of the corresponding functional at the distribution in the argument is computed.
}
\section{Methods}{
\describe{

   \item{\code{Qn}, \code{signature(x = "Any")}:}{
    interface to the \pkg{robustbase}-function \code{Qn} --- see \code{\link[robustbase]{Qn}}.}
   \item{\code{Qn}, \code{signature(x = "UnivariateDistribution")}:}{
    Qn of univariate distributions.}
   \item{\code{Qn}, \code{signature(x = "DiscreteDistribution")}:}{
    Qn of discrete distributions.}
   \item{\code{Qn}, \code{signature(x = "AffLinDistribution")}:}{\code{abs(x@a) * Qn(x@X0)}}

   \item{\code{Sn}, \code{signature(x = "Any")}:}{
    interface to the \pkg{robustbase}-function \code{Qn} --- see \code{\link[robustbase]{Sn}}.}
   \item{\code{Sn}, \code{signature(x = "UnivariateDistribution")}:}{
    Sn of univariate distributions using pseudo-random variables (Thx to N. Horbenko).}
   \item{\code{Sn}, \code{signature(x = "DiscreteDistribution")}:}{
    Sn of discrete distributions.}
   \item{\code{Sn}, \code{signature(x = "AffLinDistribution")}:}{\code{abs(x@a) * Sn(x@X0)}}

  \item{\code{var}, \code{signature(x = "Gumbel")}:}{
    exact evaluation using explicit expressions.}
  \item{\code{var}, \code{signature(x = "GPareto")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{var}, \code{signature(x = "GEV")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{var}, \code{signature(x = "Pareto")}:}{
    exact evaluation using explicit expressions.}

  \item{\code{IQR}, \code{signature(x = "Gumbel")}:}{
    exact evaluation using explicit expressions.}
  \item{\code{IQR}, \code{signature(x = "GPareto")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{IQR}, \code{signature(x = "GEV")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{IQR}, \code{signature(x = "Pareto")}:}{
    exact evaluation using explicit expressions.}

  \item{\code{median}, \code{signature(x = "Gumbel")}:}{
    exact evaluation using explicit expressions.}
  \item{\code{median}, \code{signature(x = "GEV")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{median}, \code{signature(x = "GPareto")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{median}, \code{signature(x = "Pareto")}:}{
    exact evaluation using explicit expressions.}

  \item{\code{skewness}, \code{signature(x = "Gumbel")}:}{
    exact evaluation using explicit expressions.}
  \item{\code{skewness}, \code{signature(x = "GEV")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{skewness}, \code{signature(x = "GPareto")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{skewness}, \code{signature(x = "Pareto")}:}{
    exact evaluation using explicit expressions.}

  \item{\code{kurtosis}, \code{signature(x = "Gumbel")}:}{
    exact evaluation using explicit expressions.}
  \item{\code{kurtosis}, \code{signature(x = "GEV")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{kurtosis}, \code{signature(x = "GPareto")}:}{ 
    exact evaluation using explicit expressions.}
  \item{\code{kurtosis}, \code{signature(x = "Pareto")}:}{
    exact evaluation using explicit expressions.}

  \item{\code{Sn}, \code{signature(x = "Norm")}:}{
    exact evaluation using explicit expressions.}
  \item{\code{Sn}, \code{signature(x = "GPareto")}:}{
    speeded up using interpolation grid.}
  \item{\code{Sn}, \code{signature(x = "GEV")}:}{
    speeded up using interpolation grid.}
  \item{\code{Sn}, \code{signature(x = "Gammad")}:}{
    speeded up using interpolation grid.}
  \item{\code{Sn}, \code{signature(x = "Weibull")}:}{
    speeded up using interpolation grid.}
  \item{\code{Sn}, \code{signature(x = "Pareto")}:}{
    speeded up using interpolation grid.}

  \item{\code{Qn}, \code{signature(x = "Norm")}:}{
    exact evaluation using explicit expressions.}
}}
%\references{ ~put references to the literature/web site here ~ }
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{ ~~further notes~~ }

\section{Caveat}{
If any of the packages \pkg{e1071}, \pkg{moments}, \pkg{fBasics} is to be used together with \pkg{distrEx} 
(or \pkg{RobExtremes}) the latter must be attached \emph{after} any of the first mentioned.
Otherwise \code{kurtosis()} and \code{skewness()}
defined as \emph{methods} in \pkg{distrEx} (or \pkg{RobExtremes}) may get masked.\cr
To re-mask, you may use 
\code{kurtosis <- distrEx::kurtosis; skewness <- distrEx::skewness}. 
See also \code{distrExMASK()}.
}

\examples{
# Variance of Exp(1) distribution
G <- GPareto()
var(G)

#median(Exp())
IQR(G)

## note the timing
system.time(print(Sn(GPareto(shape=0.5,scale=2))))
\donttest{
system.time(print(Sn(as(GPareto(shape=0.5,scale=2),"AbscontDistribution"))))
}
}
\seealso{\code{\link[distrEx]{Var}},\cr
 \code{\link[stats]{sd}}, \code{\link[stats:cor]{var}}, \code{\link[stats]{IQR}},\cr
 \code{\link[stats]{median}}, \code{\link[stats]{mad}},  \code{\link[distr:sd-methods]{sd}},\cr
 \code{\link[robustbase]{Sn}},  \code{\link[robustbase]{Qn}}}
\concept{functional}
\concept{var}
\concept{IQR}
\concept{median}
\concept{skewness}
\concept{kurtosis}
\concept{Sn}
\concept{Qn}
\keyword{methods}
\keyword{distribution}
\concept{integration}
