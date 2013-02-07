\name{Gumbel-class}
\docType{class}
\alias{Gumbel-class}
\alias{initialize,Gumbel-method}
\alias{loc,Gumbel-method}
\alias{loc<-,Gumbel-method}
\alias{scale,Gumbel-method}
\alias{scale<-,Gumbel-method}
\alias{+,Gumbel,numeric-method}
\alias{*,Gumbel,numeric-method}
\alias{E}
\alias{E-methods}
\alias{E,Gumbel,missing,missing-method}
\alias{var}
\alias{var-methods}
\alias{var,Gumbel-method}
\alias{skewness}
\alias{skewness-methods}
\alias{skewness,Gumbel-method}
\alias{kurtosis}
\alias{kurtosis-methods}
\alias{kurtosis,Gumbel-method}

\title{Gumbel distribution}
\description{The Gumbel cumulative distribution function with 
  location parameter \code{loc} \eqn{= \mu}{= mu} and scale 
  parameter \code{scale} \eqn{= \sigma}{= sigma} is
  \deqn{F(x) = \exp(-\exp[-(x-\mu)/\sigma])}{F(x) = exp(-exp[-(x-mu)/sigma])}
  for all real x, where \eqn{\sigma > 0}{sigma > 0}; 
  c.f. \code{rgumbel}. This distribution is also known as
  extreme value distribution of type I; confer Chapter~22 of 
  Johnson et al. (1995).
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Gumbel", loc, scale)}.
  More frequently they are created via the generating function 
  \code{Gumbel}.
}
\section{Slots}{
  \describe{
    \item{\code{img}}{Object of class \code{"Reals"}. }
    \item{\code{param}}{Object of class \code{"GumbelParameter"}. }
    \item{\code{r}}{\code{rgumbel}}
    \item{\code{d}}{\code{dgumbel}}
    \item{\code{p}}{\code{pgumbel}}
    \item{\code{q}}{\code{qgumbel}}
    \item{\code{gaps}}{(numeric) matrix or \code{NULL}}
    \item{\code{.withArith}}{logical: used internally to issue warnings as to 
            interpretation of arithmetics}
    \item{\code{.withSim}}{logical: used internally to issue warnings as to 
          accuracy}
    \item{\code{.logExact}}{logical: used internally to flag the case where 
    there are explicit formulae for the log version of density, cdf, and 
    quantile function}
    \item{\code{.lowerExact}}{logical: used internally to flag the case where 
    there are explicit formulae for the lower tail version of cdf and quantile 
    function}
    \item{\code{Symmetry}}{object of class \code{"DistributionSymmetry"};
     used internally to avoid unnecessary calculations.}
  }
}
\section{Extends}{
Class \code{"AbscontDistribution"}, directly.\cr
Class \code{"UnivariateDistribution"}, by class \code{"AbscontDistribution"}.\cr
Class \code{"Distribution"}, by class \code{"AbscontDistribution"}.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Gumbel")}: initialize method. }

    \item{loc}{\code{signature(object = "Gumbel")}: wrapped access method for
      slot \code{loc} of slot \code{param}. }

    \item{scale}{\code{signature(x = "Gumbel")}: wrapped access method for
      slot \code{scale} of slot \code{param}. }

    \item{loc<-}{\code{signature(object = "Gumbel")}: wrapped replace method for
      slot \code{loc} of slot \code{param}. }

    \item{scale<-}{\code{signature(x = "Gumbel")}: wrapped replace method for
      slot \code{scale} of slot \code{param}. }

    \item{\code{+}}{\code{signature(e1 = "Gumbel", e2 = "numeric")}: result again of 
      class \code{"Gumbel"}; exact. }

    \item{\code{*}}{\code{signature(e1 = "Gumbel", e2 = "numeric")}: result again of 
      class \code{"Gumbel"}; exact. }

    \item{E}{\code{signature(object = "Gumbel", fun = "missing", cond = "missing")}: 
      exact evaluation of expectation using explicit expressions.}

     \item{var}{\code{signature(x = "Gumbel")}: 
      exact evaluation of expectation using explicit expressions.}

     \item{skewness}{\code{signature(x = "Gumbel")}: 
      exact evaluation of expectation using explicit expressions.}

     \item{kurtosis}{\code{signature(x = "Gumbel")}: 
      exact evaluation of expectation using explicit expressions.}

     \item{median}{\code{signature(x = "Gumbel")}: 
      exact evaluation of expectation using explicit expressions.}

     \item{IQR}{\code{signature(x = "Gumbel")}: 
      exact evaluation of expectation using explicit expressions.}
  }
}
\usage{
E(object, fun, cond, ...)
\S4method{E}{Gumbel,missing,missing}(object, low = NULL, upp = NULL, ...)
var(x, ...)
\S4method{var}{Gumbel}(x, ...)
skewness(x, ...)
\S4method{skewness}{Gumbel}(x, ...)
kurtosis(x, ...)
\S4method{kurtosis}{Gumbel}(x, ...)

}
\arguments{
  \item{object}{ object of class \code{"Distribution"}}
  \item{fun}{ if missing the (conditional) expectation is computed
    else the (conditional) expection of \code{fun} is computed. }
  \item{cond}{ if not missing the conditional expectation
    given \code{cond} is computed. }
  \item{low}{lower bound of integration range.}
  \item{upp}{upper bound of integration range.}
  \item{x}{ object of class \code{"UnivariateDistribution"}}
  \item{\dots}{ additional arguments to \code{fun} }
}


\references{Johnson et al. (1995) \emph{Continuous Univariate Distributions. Vol. 2. 2nd ed.}
    New York: Wiley.}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{This class is based on the code provided by the package \pkg{evd}.}
\seealso{\code{\link[evd:gumbel]{rgumbel}}, \code{\link[distr]{AbscontDistribution-class}}}
\examples{
(G1 <- new("Gumbel", loc = 1, scale = 2))
plot(G1)
loc(G1)
scale(G1)
loc(G1) <- -1
scale(G1) <- 2
plot(G1)
}
\concept{Gumbel}
\keyword{distribution}
\concept{extreme value distribution}
\concept{absolutely continuous distribution}
\concept{S4 distribution class}
