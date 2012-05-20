\name{E}
\alias{E}
\alias{E-methods}
\alias{E,UnivariateDistribution,missing,missing-method}
\alias{E,AbscontDistribution,missing,missing-method}
\alias{E,DiscreteDistribution,missing,missing-method}
\alias{E,LatticeDistribution,missing,missing-method}
\alias{E,AffLinDistribution,missing,missing-method}
\alias{E,AffLinAbscontDistribution,missing,missing-method}
\alias{E,AffLinDiscreteDistribution,missing,missing-method}
\alias{E,AffLinLatticeDistribution,missing,missing-method}
\alias{E,MultivariateDistribution,missing,missing-method}
\alias{E,DiscreteMVDistribution,missing,missing-method}
\alias{E,UnivarLebDecDistribution,missing,missing-method}
\alias{E,AffLinUnivarLebDecDistribution,missing,missing-method}
\alias{E,UnivarMixingDistribution,missing,missing-method}
\alias{E,UnivariateDistribution,function,missing-method}
\alias{E,AbscontDistribution,function,missing-method}
\alias{E,DiscreteDistribution,function,missing-method}
\alias{E,LatticeDistribution,function,missing-method}
\alias{E,MultivariateDistribution,function,missing-method}
\alias{E,DiscreteMVDistribution,function,missing-method}
\alias{E,UnivarLebDecDistribution,function,missing-method}
\alias{E,UnivarMixingDistribution,function,missing-method}
\alias{E,AcDcLcDistribution,ANY,ANY-method}
\alias{E,CompoundDistribution,missing,missing-method}
\alias{E,UnivariateCondDistribution,missing,numeric-method}
\alias{E,AbscontCondDistribution,missing,numeric-method}
\alias{E,DiscreteCondDistribution,missing,numeric-method}
\alias{E,UnivarLebDecDistribution,missing,ANY-method}
\alias{E,UnivarMixingDistribution,missing,ANY-method}
\alias{E,UnivarLebDecDistribution,function,ANY-method}
\alias{E,UnivariateCondDistribution,function,numeric-method}
\alias{E,UnivarMixingDistribution,function,ANY-method}
\alias{E,AbscontCondDistribution,function,numeric-method}
\alias{E,DiscreteCondDistribution,function,numeric-method}
\alias{E,Arcsine,missing,missing-method}
\alias{E,Beta,missing,missing-method}
\alias{E,Binom,missing,missing-method}
\alias{E,Cauchy,missing,missing-method}
\alias{E,Chisq,missing,missing-method}
\alias{E,Dirac,missing,missing-method}
\alias{E,DExp,missing,missing-method}
\alias{E,Exp,missing,missing-method}
\alias{E,Fd,missing,missing-method}
\alias{E,Gammad,missing,missing-method}
\alias{E,Gammad,function,missing-method}
\alias{E,Geom,missing,missing-method}
\alias{E,Gumbel,missing,missing-method}
\alias{E,GPareto,missing,missing-method}
\alias{E,GPareto,function,missing-method}
\alias{E,GEV,missing,missing-method}
\alias{E,GEV,function,missing-method}
\alias{E,Hyper,missing,missing-method}
\alias{E,Logis,missing,missing-method}
\alias{E,Lnorm,missing,missing-method}
\alias{E,Nbinom,missing,missing-method}
\alias{E,Norm,missing,missing-method}
\alias{E,Pareto,missing,missing-method}
\alias{E,Pois,missing,missing-method}
\alias{E,Td,missing,missing-method}
\alias{E,Unif,missing,missing-method}
\alias{E,Weibull,missing,missing-method}

\title{Generic Function for the Computation of (Conditional) Expectations}
\description{
  Generic function for the computation of (conditional) expectations.
}
\usage{
E(object, fun, cond, ...)

\S4method{E}{UnivariateDistribution,missing,missing}(object, 
             low = NULL, upp = NULL, Nsim = getdistrExOption("MCIterations"), ...)

\S4method{E}{UnivariateDistribution,function,missing}(object, fun, 
        useApply = TRUE, low = NULL, upp = NULL,  Nsim = getdistrExOption("MCIterations"), ...)

\S4method{E}{AbscontDistribution,function,missing}(object, fun, useApply = TRUE,
             low = NULL, upp = NULL, 
             rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ...)

\S4method{E}{UnivarMixingDistribution,missing,missing}(object, low = NULL, 
             upp = NULL, rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ...)

\S4method{E}{UnivarMixingDistribution,function,missing}(object, fun, useApply = TRUE, 
             low = NULL, upp = NULL, 
             rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ...)

\S4method{E}{UnivarMixingDistribution,missing,ANY}(object, cond, low = NULL, 
             upp = NULL, rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ...)

\S4method{E}{UnivarMixingDistribution,function,ANY}(object, fun, cond, useApply = TRUE,
             low = NULL, upp = NULL, rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ...)

\S4method{E}{DiscreteDistribution,function,missing}(object, fun, useApply = TRUE, 
             low = NULL, upp = NULL, ...)

\S4method{E}{AffLinDistribution,missing,missing}(object, low = NULL, upp = NULL, ...)

\S4method{E}{AffLinUnivarLebDecDistribution,missing,missing}(object, low = NULL, upp = NULL, ...)

\S4method{E}{MultivariateDistribution,missing,missing}(object, 
             Nsim = getdistrExOption("MCIterations"), ...)
\S4method{E}{MultivariateDistribution,function,missing}(object, fun, useApply = TRUE, 
             Nsim = getdistrExOption("MCIterations"), ...)

\S4method{E}{DiscreteMVDistribution,missing,missing}(object, low = NULL, upp = NULL, ...)

\S4method{E}{DiscreteMVDistribution,function,missing}(object, fun, 
             useApply = TRUE, ...)

\S4method{E}{AbscontCondDistribution,missing,numeric}(object, cond, useApply = TRUE,
             low = NULL, upp = NULL, 
             rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ...)

\S4method{E}{DiscreteCondDistribution,missing,numeric}(object, cond, useApply = TRUE,
             low = NULL, upp = NULL, ...)

\S4method{E}{UnivariateCondDistribution,function,numeric}(object, fun, cond, 
              withCond = FALSE, useApply = TRUE, low = NULL, upp = NULL,
              Nsim = getdistrExOption("MCIterations"), ...)

\S4method{E}{AbscontCondDistribution,function,numeric}(object, fun, cond, 
               withCond = FALSE, useApply = TRUE, low = NULL, upp = NULL,
             rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac")
             , ...)

\S4method{E}{DiscreteCondDistribution,function,numeric}(object, fun, cond, 
             withCond = FALSE, useApply = TRUE, low = NULL, upp = NULL,...)

\S4method{E}{DiscreteCondDistribution,function,numeric}(object, fun, cond, 
             withCond = FALSE, useApply = TRUE, low = NULL, upp = NULL,...)
             
\S4method{E}{UnivarLebDecDistribution,missing,missing}(object, low = NULL, upp = NULL,  
             rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ... )
\S4method{E}{UnivarLebDecDistribution,function,missing}(object, fun, 
             useApply = TRUE, low = NULL, upp = NULL, rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ... )
\S4method{E}{UnivarLebDecDistribution,missing,ANY}(object, cond, 
             low = NULL, upp = NULL, rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ... )
\S4method{E}{UnivarLebDecDistribution,function,ANY}(object, fun, cond, 
             useApply = TRUE, low = NULL, upp = NULL,
             rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ... )

\S4method{E}{AcDcLcDistribution,ANY,ANY}(object, fun, cond, 
             low = NULL, upp = NULL, rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac"), ... )
\S4method{E}{CompoundDistribution,missing,missing}(object, low = NULL, upp = NULL, ...)

\S4method{E}{Arcsine,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Beta,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Binom,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Cauchy,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Chisq,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Dirac,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{DExp,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Exp,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Fd,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Gammad,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Gammad,function,missing}(object, fun, low = NULL, upp = NULL, 
                                      rel.tol = getdistrExOption("ErelativeTolerance"),
                                      lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
                                      upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
                                      IQR.fac = getdistrExOption("IQR.fac"), ...)
\S4method{E}{Geom,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Gumbel,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{GPareto,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{GPareto,function,missing}(object, fun, low = NULL, upp = NULL, 
                                      rel.tol = getdistrExOption("ErelativeTolerance"),
                                      lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"),
                                      upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
                                      IQR.fac = max(10000, getdistrExOption("IQR.fac")), ...)
\S4method{E}{Hyper,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Logis,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Lnorm,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Nbinom,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Norm,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Pareto,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Pois,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Unif,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Td,missing,missing}(object, low = NULL, upp = NULL, ...)
\S4method{E}{Weibull,missing,missing}(object, low = NULL, upp = NULL, ...)

}
\arguments{
  \item{object}{ object of class \code{"Distribution"}}
  \item{fun}{ if missing the (conditional) expectation is computed
    else the (conditional) expection of \code{fun} is computed. }
  \item{cond}{ if not missing the conditional expectation 
    given \code{cond} is computed. }
  \item{Nsim}{ number of MC simulations used to determine the expectation. }
  \item{rel.tol}{relative tolerance for \code{distrExIntegrate}.} 
  \item{low}{lower bound of integration range.}
  \item{upp}{upper bound of integration range.}
  \item{lowerTruncQuantile}{lower quantile for quantile based integration range.}
  \item{upperTruncQuantile}{upper quantile for quantile based integration range.}
  \item{IQR.fac}{factor for scale based integration range (i.e.; 
  median of the distribution \eqn{\pm}{+-}\code{IQR.fac}\eqn{\times}{*}IQR).}
  \item{\dots}{ additional arguments to \code{fun} }
  \item{useApply}{ logical: should \code{sapply}, respectively \code{apply} 
    be used to evaluate \code{fun}. }
  \item{withCond}{ logical: is \code{cond} in the argument list of \code{fun}. }
}
\details{The precision of the computations can be controlled via 
  certain global options; cf. \code{\link{distrExOptions}}. 
    Also note that arguments \code{low} and \code{upp} should be given as
  named arguments in order to prevent them to be matched by arguments
  \code{fun} or \code{cond}. Also the result, when arguments 
  \code{low} or \code{upp} is given, is the \emph{unconditional value} of the
  expectation; no conditioning with respect to \code{low <= object <= upp}
  is done.}

\value{
  The (conditional) expectation is computed.
}
\section{Methods}{
\describe{
  \item{object = "UnivariateDistribution", fun = "missing", cond = "missing":}{ 
    expectation of univariate distributions using crude Monte-Carlo integration. }

  \item{object = "AbscontDistribution", fun = "missing", cond = "missing":}{ 
    expectation of absolutely continuous univariate distributions
    using \code{distrExIntegrate}. }

  \item{object = "DiscreteDistribution", fun = "missing", cond = "missing":}{  
    expectation of discrete univariate distributions using \code{support}
    and \code{sum}.}

  \item{object = "MultivariateDistribution", fun = "missing", cond = "missing":}{
    expectation of multivariate distributions using crude Monte-Carlo integration. }

  \item{object = "DiscreteMVDistribution", fun = "missing", cond = "missing":}{
    expectation of discrete multivariate distributions. The computation is based
    on \code{support} and \code{sum}.}

  \item{object = "UnivariateDistribution", fun = "missing", cond = "missing":}{ 
    expectation of univariate Lebesgue decomposed distributions
    by separate calculations for discrete and absolutely continuous part. }

  \item{object = "AffLinDistribution", fun = "missing", cond = "missing":}{
        expectation of an affine linear transformation \eqn{aX+b} as
        \eqn{a E[X]+b} for \code{X} either \code{"DiscreteDistribution"}
        or \code{"AbscontDistribution"}.
        }
  \item{object = "AffLinUnivarLebDecDistribution", fun = "missing", cond = "missing":}{
        expectation of an affine linear transformation \eqn{aX+b} as
        \eqn{a E[X]+b} for \code{X} either \code{"UnivarLebDecDistribution"}.
        }

  \item{object = "UnivariateDistribution", fun = "function", cond = "missing":}{ 
    expectation of \code{fun} under univariate distributions using 
    crude Monte-Carlo integration. }

  \item{object = "UnivariateDistribution", fun = "function", cond = "missing":}{ 
    expectation of \code{fun} under univariate Lebesgue decomposed distributions
    by separate calculations for discrete and absolutely continuous part. }

  \item{object = "AbscontDistribution", fun = "function", cond = "missing":}{ 
    expectation of \code{fun} under absolutely continuous 
    univariate distributions using \code{distrExIntegrate}. }

  \item{object = "DiscreteDistribution", fun = "function", cond = "missing":}{  
    expectation of \code{fun} under discrete univariate 
    distributions using \code{support} and \code{sum}. }

  \item{object = "MultivariateDistribution", fun = "function", cond = "missing":}{
    expectation of multivariate distributions using crude Monte-Carlo integration. }

  \item{object = "DiscreteMVDistribution", fun = "function", cond = "missing":}{
    expectation of \code{fun} under discrete multivariate 
    distributions. The computation is based on \code{support} and \code{sum}. }  

  \item{object = "UnivariateCondDistribution", fun = "missing", cond = "numeric":}{
    conditional expectation for univariate conditional distributions given \code{cond}. 
    The integral is computed using crude Monte-Carlo integration. }

  \item{object = "AbscontCondDistribution", fun = "missing", cond = "numeric":}{
    conditional expectation for absolutely continuous, univariate 
    conditional distributions given \code{cond}. The computation
    is based on \code{distrExIntegrate}. }

  \item{object = "DiscreteCondDistribution", fun = "missing", cond = "numeric":}{ 
    conditional expectation for discrete, univariate conditional 
    distributions given \code{cond}. The computation is based 
    on \code{support} and \code{sum}. }

  \item{object = "UnivariateCondDistribution", fun = "function", cond = "numeric":}{
    conditional expectation of \code{fun} under univariate conditional distributions 
    given \code{cond}. The integral is computed using crude Monte-Carlo integration. }

  \item{object = "AbscontCondDistribution", fun = "function", cond = "numeric":}{
    conditional expectation of \code{fun} under absolutely continuous, 
    univariate conditional distributions given \code{cond}. The
    computation is based on \code{distrExIntegrate}. }

  \item{object = "DiscreteCondDistribution", fun = "function", cond = "numeric":}{ 
    conditional expectation of \code{fun} under discrete, univariate 
    conditional distributions given \code{cond}. The computation is
    based on \code{support} and \code{sum}. }

  \item{object = "UnivarLebDecDistribution", fun = "missing", cond = "missing":}{
        expectation by separate evaluation of expectation of discrete and
        abs. continuous part and subsequent weighting.
        }
  \item{object = "UnivarLebDecDistribution", fun = "function", cond = "missing":}{
        expectation by separate evaluation of expectation of discrete and
        abs. continuous part and subsequent weighting.
        }
  \item{object = "UnivarLebDecDistribution", fun = "missing", cond = "ANY":}{
        expectation by separate evaluation of expectation of discrete and
        abs. continuous part and subsequent weighting.
        }
  \item{object = "UnivarLebDecDistribution", fun = "function", cond = "ANY":}{
        expectation by separate evaluation of expectation of discrete and
        abs. continuous part and subsequent weighting.
        }
  \item{object = "UnivarMixingDistribution", fun = "missing", cond = "missing":}{
        expectation is computed component-wise with subsequent weighting acc.
        to \code{mixCoeff}.}
  \item{object = "UnivarMixingDistribution", fun = "function", cond = "missing":}{
        expectation is computed component-wise with subsequent weighting acc.
        to \code{mixCoeff}.}
  \item{object = "UnivarMixingDistribution", fun = "missing", cond = "ANY":}{
        expectation is computed component-wise with subsequent weighting acc.
        to \code{mixCoeff}.}
  \item{object = "UnivarMixingDistribution", fun = "function", cond = "ANY":}{
        expectation is computed component-wise with subsequent weighting acc.
        to \code{mixCoeff}.}

  \item{object = "AcDcLcDistribution", fun = "ANY", cond = "ANY":}{
        expectation by first coercing to class \code{"UnivarLebDecDistribution"}
        and using the corresponding method.
        }
  \item{object = "CompoundDistribution", fun = "missing", cond = "missing":}{
        if we are in i.i.d. situation (i.e., slot \code{SummandsDistr} is of
        class \code{UnivariateDistribution}) the formula 
        \eqn{E[N]E[S]} for \eqn{N} the frequency distribution and
        \eqn{S} the summand distribution; else we coerce to 
        \code{"UnivarLebDecDistribution"}.
        }

  \item{object = "Arcsine", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Beta", fun = "missing", cond = "missing":}{ 
    for noncentrality 0 exact evaluation using explicit expressions.}
  \item{object = "Binom", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Cauchy", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Chisq", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Dirac", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "DExp", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Exp", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Fd", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Gammad", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Gammad", fun = "function", cond = "missing":}{ 
    use substitution method (y := log(x)) for numerical integration.}
  \item{object = "Geom", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Gumbel", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "GPareto", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "GPareto", fun = "function", cond = "missing":}{ 
    use substitution method (y := log(x)) for numerical integration.}
  \item{object = "Hyper", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Logis", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Lnorm", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Nbinom", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Norm", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Pareto", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Pois", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Unif", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Td", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
  \item{object = "Weibull", fun = "missing", cond = "missing":}{ 
    exact evaluation using explicit expressions.}
}}
%\references{ ~put references to the literature/web site here ~ }
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de} and Peter Ruckdeschel \email{peter.ruckdeschel@uni-bayreuth.de}}
%\note{ ~~further notes~~ }
\seealso{\code{\link{distrExIntegrate}}, \code{\link{m1df}}, \code{\link{m2df}},
    \code{\link[distr]{Distribution-class}}}
\examples{
# mean of Exp(1) distribution
E <- Exp() 

E(E) ## uses explicit terms
E(as(E,"AbscontDistribution")) ## uses numerical integration
E(as(E,"UnivariateDistribution")) ## uses simulations
E(E, fun = function(x){2*x^2}) ## uses simulations

# the same operator for discrete distributions:
P <- Pois(lambda=2)

E(P) ## uses explicit terms
E(as(P,"DiscreteDistribution")) ## uses sums
E(as(P,"UnivariateDistribution")) ## uses simulations
E(P, fun = function(x){2*x^2}) ## uses simulations


# second moment of N(1,4)
E(Norm(mean=1, sd=2), fun = function(x){x^2})
E(Norm(mean=1, sd=2), fun = function(x){x^2}, useApply = FALSE)

# conditional distribution of a linear model
D1 <- LMCondDistribution(theta = 1) 
E(D1, cond = 1)
E(Norm(mean=1))
E(D1, function(x){x^2}, cond = 1)
E(Norm(mean=1), fun = function(x){x^2})
E(D1, function(x, cond){cond*x^2}, cond = 2, withCond = TRUE, useApply = FALSE)
E(Norm(mean=2), function(x){2*x^2})

E(as(Norm(mean=2),"AbscontDistribution"))
### somewhat less accurate:
E(as(Norm(mean=2),"AbscontDistribution"), 
     lowerTruncQuantil=1e-4,upperTruncQuantil=1e-4, IQR.fac= 4)
### even less accurate:
E(as(Norm(mean=2),"AbscontDistribution"), 
     lowerTruncQuantil=1e-2,upperTruncQuantil=1e-2, IQR.fac= 4)
### no good idea, but just as an example:
E(as(Norm(mean=2),"AbscontDistribution"), 
     lowerTruncQuantil=1e-2,upperTruncQuantil=1e-2, IQR.fac= .1)

### truncation of integration range; see also m1df...
E(Norm(mean=2), low=2,upp=4)

E(Cauchy())
E(Cauchy(),upp=3,low=-2)
# some Lebesgue decomposed distribution 
mymix <- UnivarLebDecDistribution(acPart = Norm(), discretePart = Binom(4,.4),
         acWeight = 0.4)
E(mymix)
}
\concept{expectation}
\concept{conditional expectation}
\concept{integration}
\keyword{methods}
\keyword{distribution}
\concept{functional}
\concept{moment}
