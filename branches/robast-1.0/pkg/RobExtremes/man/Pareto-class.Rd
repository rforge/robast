\name{Pareto-class}
\docType{class}
\alias{Pareto-class}
\alias{initialize,Pareto-method}
\alias{shape,Pareto-method}
\alias{shape<-,Pareto-method}
\alias{Min,Pareto-method}
\alias{Min<-,Pareto-method}

\title{Pareto distribution}
\description{[borrowed from \pkg{actuar}]:
 
  The (Single-parameter) Pareto distribution with parameter \code{shape}
  \eqn{= \alpha}{= a} has density:
  \deqn{f(x) = \frac{\alpha \theta^\alpha}{x^{\alpha + 1}}}{%
    f(x) = a b^a/x^(a + 1)}
  for \eqn{x > \theta}{x > b}, \eqn{\alpha > 0}{a > 0} and \eqn{\theta >
    0}{b > 0}.
  
  Although there appears to be two parameters, only \code{shape} is a true
  parameter. The value of \code{min} \eqn{= \theta}{= b} must be set in
  advance.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Pareto", shape, Min)}.
  More frequently they are created via the generating function 
  \code{Pareto}.
}
\section{Slots}{
  \describe{
    \item{\code{img}}{Object of class \code{"Reals"}. }
    \item{\code{param}}{Object of class \code{"ParetoParameter"}. }
    \item{\code{r}}{\code{rpareto1}}
    \item{\code{d}}{\code{dpareto1}}
    \item{\code{p}}{\code{ppareto1}}
    \item{\code{q}}{\code{qpareto1}}
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
    \item{initialize}{\code{signature(.Object = "Pareto")}: initialize method. }

    \item{shape}{\code{signature(object = "Pareto")}: wrapped access method for
      slot \code{shape} of slot \code{param}. }

    \item{Min}{\code{signature(x = "Pareto")}: wrapped access method for
      slot \code{Min} of slot \code{param}. }

    \item{shape<-}{\code{signature(object = "Pareto")}: wrapped replace method for
      slot \code{shape} of slot \code{param}. }

    \item{Min<-}{\code{signature(x = "Pareto")}: wrapped replace method for
      slot \code{Min} of slot \code{param}. }

    \item{E}{\code{signature(object = "Pareto", fun = "missing", cond = "missing")}: 
    exact evaluation using explicit expressions.}
    \item{var}{\code{signature(signature(x = "Pareto")}: 
    exact evaluation using explicit expressions.}
    \item{median}{\code{signature(signature(x = "Pareto")}: 
    exact evaluation using explicit expressions.}
    \item{IQR}{\code{signature(signature(x = "Pareto")}:
    exact evaluation using explicit expressions.}
    \item{skewness}{\code{signature(signature(x = "Pareto")}: 
    exact evaluation using explicit expressions.}
    \item{kurtosis}{\code{signature(signature(x = "Pareto")}: 
    exact evaluation using explicit expressions.}
  }
}
\references{Johnson et al. (1995) \emph{Continuous Univariate Distributions. Vol. 2. 2nd ed.}
    New York: Wiley.\cr
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
    }
\author{Nataliya Horbenko \email{Nataliya.Horbenko@itwm.fraunhofer.de}}
\note{This class is based on the code provided by the package \pkg{actuar}
by  Vincent Goulet and Mathieu Pigeon.}
\seealso{\code{\link[actuar:SingleParameterPareto]{dpareto1}}, \code{\link[distr]{AbscontDistribution-class}}}
\examples{
(P1 <- new("Pareto", shape = 1, Min = 2))
plot(P1)
shape(P1)
Min(P1)
shape(P1) <- 4
Min(P1) <- 2
plot(P1)
}
\concept{Pareto}
\keyword{distribution}
\concept{extreme value distribution}
\concept{absolutely continuous distribution}
\concept{S4 distribution class}
