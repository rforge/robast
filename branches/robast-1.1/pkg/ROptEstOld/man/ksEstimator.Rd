\name{ksEstimator}
\alias{ksEstimator}
\alias{ksEstimator-methods}
\alias{ksEstimator,numeric,Binom-method}
\alias{ksEstimator,numeric,Pois-method}
\alias{ksEstimator,numeric,Norm-method}
\alias{ksEstimator,numeric,Lnorm-method}
\alias{ksEstimator,numeric,Gumbel-method}
\alias{ksEstimator,numeric,Exp-method}
\alias{ksEstimator,numeric,Gammad-method}

\title{Generic Function for the Computation of the Kolmogorov Minimum Distance Estimator}
\description{
  Generic function for the computation of the Kolmogorov(-Smirnov)
  minimum distance estimator.
}
\usage{
ksEstimator(x, distribution, ...)

\S4method{ksEstimator}{numeric,Binom}(x, distribution, param, eps = .Machine$double.eps^0.5)

\S4method{ksEstimator}{numeric,Pois}(x, distribution, param, eps = .Machine$double.eps^0.5)

\S4method{ksEstimator}{numeric,Norm}(x, distribution, param, eps = .Machine$double.eps^0.5)

\S4method{ksEstimator}{numeric,Lnorm}(x, distribution, param, eps = .Machine$double.eps^0.5)

\S4method{ksEstimator}{numeric,Gumbel}(x, distribution, param, eps = .Machine$double.eps^0.5)

\S4method{ksEstimator}{numeric,Exp}(x, distribution, param, eps = .Machine$double.eps^0.5)

\S4method{ksEstimator}{numeric,Gammad}(x, distribution, param, eps = .Machine$double.eps^0.5)
}
\arguments{
  \item{x}{ sample }
  \item{distribution}{ object of class \code{"Distribution"} }
  \item{\dots}{ additional parameters }
  \item{param}{ name of the unknown parameter. If missing all parameters
    of the corresponding distribution are estimated. }
  \item{eps}{ the desired accuracy (convergence tolerance). }
}
\details{In case of discrete distributions the Kolmogorov distance is computed and
  the parameters which lead to the minimum distance are returned. In case of 
  absolutely continuous distributions \code{ks.test} is called and the parameters
  which minimize the corresponding test statistic are returned. }
\value{The Kolmogorov minimum distance estimator is computed. Returns a list
  with components named like the parameters of \code{distribution}. }
\section{Methods}{
\describe{
  \item{x = "numeric", distribution = "Binom"}{ Binomial distributions. }

  \item{x = "numeric", distribution = "Pois"}{ Poisson distributions. }

  \item{x = "numeric", distribution = "Norm"}{ Normal distributions. }

  \item{x = "numeric", distribution = "Lnorm"}{ Lognormal distributions. }

  \item{x = "numeric", distribution = "Gumbel"}{ Gumbel distributions. }

  \item{x = "numeric", distribution = "Exp"}{ Exponential distributions. }

  \item{x = "numeric", distribution = "Gamma"}{ Gamma distributions. }
}}
\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distr]{Distribution-class}}}
\examples{
x <- rnorm(100, mean = 1, sd = 2)
ksEstimator(x=x, distribution = Norm()) # estimate mean and sd
ksEstimator(x=x, distribution = Norm(mean = 1), param = "sd") # estimate sd
ksEstimator(x=x, distribution = Norm(sd = 2), param = "mean") # estimate mean
mean(x)
median(x)
sd(x)
mad(x)
}
\concept{Kolmogorov minimum distance estimator}
\concept{minimum distance estimator}
\concept{estimator}
\keyword{robust}
