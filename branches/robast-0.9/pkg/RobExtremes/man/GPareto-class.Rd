\name{GPareto-class}
\docType{class}
\alias{GPareto-class}
\alias{initialize,GPareto-method}
\alias{loc,GPareto-method}
\alias{loc<-,GPareto-method}
\alias{location,GPareto-method}
\alias{location<-,GPareto-method}
\alias{scale,GPareto-method}
\alias{scale<-,GPareto-method}
\alias{shape,GPareto-method}
\alias{shape<-,GPareto-method}
\alias{+,GPareto,numeric-method}
\alias{*,GPareto,numeric-method}


\title{Generalized Pareto distribution}
\description{[borrowed from \pkg{evd}]:
 
  The (Three-parameter) generalized Pareto distribution with parameter \code{loc}\eqn{= a},
  \code{scale} \eqn{= b}, \code{shape} \eqn{= c} has density:
  \deqn{f(x) = \frac{1}{b} (1+c z)^(-1/c - 1), \quad z = \frac{x-a}{c}}
  for \eqn{x > a} (\eqn{ c \geq 0}) and \eqn{a \leq x \leq a - b/c}(\eqn{c < 0}).

  }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("GPareto", loc, scale,shape)}.
  More frequently they are created via the generating function 
  \code{GPareto}.
}
\section{Slots}{
  \describe{
    \item{\code{img}}{Object of class \code{"Reals"}. }
    \item{\code{param}}{Object of class \code{"GParetoParameter"}. }
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
    \item{initialize}{\code{signature(.Object = "GPareto")}: initialize method. }

    \item{shape}{\code{signature(object = "GPareto")}: wrapped access method for
      slot \code{shape} of slot \code{param}. }
    \item{loc}{\code{signature(object = "GPareto")}: wrapped access method for
      slot \code{loc} of slot \code{param}. }
    \item{location}{\code{signature(object = "GPareto")}: alias to \code{loc},
      to support argument naming of package \pkg{VGAM}. }
    \item{scale}{\code{signature(x = "GPareto")}: wrapped access method for
      slot \code{scale} of slot \code{param}. }

    \item{shape<-}{\code{signature(object = "GPareto")}: wrapped replace method for
      slot \code{shape} of slot \code{param}. }
    \item{loc<-}{\code{signature(object = "GPareto")}: wrapped replace  method for
      slot \code{loc} of slot \code{param}. }
    \item{location<-}{\code{signature(object = "GPareto")}: alias to \code{loc<-},
      to support argument naming of package \pkg{VGAM}. }
    \item{scale<-}{\code{signature(x = "GPareto")}: wrapped replace  method for
      slot \code{scale} of slot \code{param}. }

    \item{+}{\code{signature(e1 = "GPareto", e2 = "numeric")}: exact method
      for this transformation --- stays within this class. }
    \item{*}{\code{signature(e1 = "GPareto", e2 = "numeric")}: exact method
      for this transformation --- stays within this class if \code{e2>0}. }

    \item{E}{\code{signature(object = "GPareto", fun = "missing", cond = "missing")}:
    exact evaluation using explicit expressions.}
    \item{var}{\code{signature(signature(x = "GPareto")}:
    exact evaluation using explicit expressions.}
    \item{median}{\code{signature(signature(x = "GPareto")}: 
    exact evaluation using explicit expressions.}
    \item{IQR}{\code{signature(signature(x = "GPareto")}: 
    exact evaluation using explicit expressions.}
    \item{skewness}{\code{signature(signature(x = "GPareto")}: 
    exact evaluation using explicit expressions.}
    \item{kurtosis}{\code{signature(signature(x = "GPareto")}:
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
(P1 <- new("GPareto", loc = 0, scale = 1,shape = 0))
plot(P1)
shape(P1)
loc(P1)
scale(P1) <- 4
loc(P1) <- 2
plot(P1)
}
\concept{GPareto}
\keyword{distribution}
\concept{extreme value distribution}
\concept{absolutely continuous distribution}
\concept{S4 distribution class}
