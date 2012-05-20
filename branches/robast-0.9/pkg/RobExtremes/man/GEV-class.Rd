\name{GEV-class}
\docType{class}
\alias{GEV-class}
\alias{initialize,GEV-method}
\alias{loc,GEV-method}
\alias{loc<-,GEV-method}
\alias{location,GEV-method}
\alias{location<-,GEV-method}
\alias{scale,GEV-method}
\alias{scale<-,GEV-method}
\alias{shape,GEV-method}
\alias{shape<-,GEV-method}
\alias{+,GEV,numeric-method}
\alias{*,GEV,numeric-method}


\title{Generalized EV distribution}
\description{[borrowed from \pkg{evd}]:
The GEV distribution function with parameters \code{loc}\eqn{= a},
  \code{scale} \eqn{= b}, \code{shape} \eqn{= s} is

\deqn{G(x) = exp[-{1+s(z-a)/b}^(-1/s)]}

for \eqn{1+s(z-a)/b > 0}, where \eqn{b > 0}. If \eqn{s = 0} the distribution is 
defined by continuity. If \eqn{1+s(z-a)/b <= 0}, the value \eqn{z} is either 
greater than the upper end point (if \eqn{s < 0}), or less than the lower end 
point (if eqn{s > 0}).
  }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("GEV", loc, scale,shape)}.
  More frequently they are created via the generating function 
  \code{GEV}.
}
\section{Slots}{
  \describe{
    \item{\code{img}}{Object of class \code{"Reals"}. }
    \item{\code{param}}{Object of class \code{"GEVParameter"}. }
    \item{\code{r}}{\code{rgpd}}
    \item{\code{d}}{\code{dgpd}}
    \item{\code{p}}{\code{pgpd}, but vectorized and with special treatment of 
             arguments \code{lower.tail} and \code{log.p}}
    \item{\code{q}}{\code{qgpd}, but vectorized and with special treatment of 
             arguments \code{lower.tail} and \code{log.p}}
    \item{\code{gaps}}{(numeric) matrix or \code{NULL}}
    \item{\code{.withArith}}{logical: used internally to issue warnings as to interpretation of arithmetics}
    \item{\code{.withSim}}{logical: used internally to issue warnings as to accuracy}
    \item{\code{.logExact}}{logical: used internally to flag the case where there are explicit formulae for the
                              log version of density, cdf, and quantile function}
    \item{\code{.lowerExact}}{logical: used internally to flag the case where there are explicit formulae for the
                              lower tail version of cdf and quantile function}
  }
}
\section{Extends}{
Class \code{"AbscontDistribution"}, directly.\cr
Class \code{"UnivariateDistribution"}, by class \code{"AbscontDistribution"}.\cr
Class \code{"Distribution"}, by class \code{"AbscontDistribution"}.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "GEV")}: initialize method. }

    \item{shape}{\code{signature(object = "GEV")}: wrapped access method for
      slot \code{shape} of slot \code{param}. }
    \item{loc}{\code{signature(object = "GEV")}: wrapped access method for
      slot \code{loc} of slot \code{param}. }
    \item{location}{\code{signature(object = "GEV")}: alias to \code{loc},
      to support argument naming of package \pkg{VGAM}. }
    \item{scale}{\code{signature(x = "GEV")}: wrapped access method for
      slot \code{scale} of slot \code{param}. }

    \item{shape<-}{\code{signature(object = "GEV")}: wrapped replace method for
      slot \code{shape} of slot \code{param}. }
    \item{loc<-}{\code{signature(object = "GEV")}: wrapped replace  method for
      slot \code{loc} of slot \code{param}. }
    \item{location<-}{\code{signature(object = "GEV")}: alias to \code{loc<-},
      to support argument naming of package \pkg{VGAM}. }
    \item{scale<-}{\code{signature(x = "GEV")}: wrapped replace  method for
      slot \code{scale} of slot \code{param}. }

    \item{+}{\code{signature(e1 = "GEV", e2 = "numeric")}: exact method
      for this transformation --- stays within this class. }
    \item{*}{\code{signature(e1 = "GEV", e2 = "numeric")}: exact method
      for this transformation --- stays within this class if \code{e2>0}. }

    \item{E}{\code{signature(object = "GEV", fun = "missing", cond = "missing")}:
    exact evaluation using explicit expressions.}
    \item{var}{\code{signature(signature(x = "GEV")}:
    exact evaluation using explicit expressions.}
    \item{median}{\code{signature(signature(x = "GEV")}: 
    exact evaluation using explicit expressions.}
    \item{IQR}{\code{signature(signature(x = "GEV")}: 
    exact evaluation using explicit expressions.}
    \item{skewness}{\code{signature(signature(x = "GEV")}: 
    exact evaluation using explicit expressions.}
    \item{kurtosis}{\code{signature(signature(x = "GEV")}:
    exact evaluation using explicit expressions.}


  }
}
\references{Pickands, J. (1975) \emph{Statistical inference using extreme order
     statistics. _Annals of Statistics_, *3*, 119-131.}}
\author{Nataliya Horbenko \email{Nataliya.Horbenko@itwm.fraunhofer.de}}
\note{This class is based on the code provided by the package \pkg{evd}
by  A. G. Stephenson.}
\seealso{\code{\link[evd:gpd]{dgpd}}, \code{\link[distr]{AbscontDistribution-class}}}
\examples{
(P1 <- new("GEV", loc = 0, scale = 1,shape = 0))
plot(P1)
shape(P1)
loc(P1)
scale(P1) <- 4
loc(P1) <- 2
plot(P1)
}
\concept{GEV}
\keyword{distribution}
\concept{extreme value distribution}
\concept{absolutely continuous distribution}
\concept{S4 distribution class}
