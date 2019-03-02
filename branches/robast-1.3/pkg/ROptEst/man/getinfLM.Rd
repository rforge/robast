\name{getInfLM}
\alias{getInfLM}
\alias{getLagrangeMultByOptim}
\alias{getLagrangeMultByIter}

\title{Functions to determine Lagrange multipliers}
\description{
  Functions to determine Lagrange multipliers \code{A} and \code{a}
  in a Hampel problem or in a(n) (inner) loop in a MSE problem; can be done
  either by optimization or by fixed point iteration. These functions are
  rarely called directly.
}
\usage{
getLagrangeMultByIter(b, L2deriv, risk, trafo,
                      neighbor, biastype, normtype, Distr,
                      a.start, z.start, A.start, w.start, std, z.comp,
                      A.comp, maxiter, tol, verbose = NULL,
                      warnit = TRUE, ...)
getLagrangeMultByOptim(b, L2deriv, risk, FI, trafo,
                      neighbor, biastype, normtype, Distr,
                      a.start, z.start, A.start, w.start,  std, z.comp,
                      A.comp, maxiter, tol, verbose = NULL, ...)

}
\arguments{
  \item{b}{numeric; (\eqn{>b_{\rm\scriptstyle min}}{>b_min}; clipping bound
  for which the Lagrange multipliers are searched}
  \item{L2deriv}{ L2-derivative of some L2-differentiable family
    of probability measures. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{FI}{ matrix: Fisher information. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{object of class \code{"BiasType"} --- the bias type with we work.}
  \item{normtype}{object of class \code{"NormType"} --- the norm type with we work.}
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{a.start}{ initial value for the centering constant (in \code{p}-space). }
  \item{z.start}{ initial value for the centering constant (in \code{k}-space). }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{w.start}{ initial value for the weight function. }
  \item{std}{ matrix of (or which may coerced to) class
                  \code{PosSemDefSymmMatrix} for use of different
                  (standardizing) norm. }
  \item{z.comp}{ logical vector: indication which components of the
    centering constant have to be computed. }
  \item{A.comp}{ matrix: indication which components of the standardizing
    matrix have to be computed. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed. }
  \item{warnit}{ logical: if \code{TRUE} warning is issued if
                  maximal number of iterations is reached. }
  \item{\dots}{ additional parameters for \code{optim} and \code{E}. }
}
%\details{}
\value{a list with items
      \item{A}{Lagrange multiplier \code{A} (standardizing matrix)}
      \item{a}{Lagrange multiplier \code{a} (centering in \code{p}-space)}
      \item{z}{Lagrange multiplier \code{z} (centering in \code{k}-space)}
      \item{w}{weight function involving Lagrange multipliers}
      \item{biastype}{(possibly modified) bias type \code{biastype} from argument}
      \item{normtype}{(possibly modified) norm type \code{normtype} from argument}
      \item{normtype.old}{(possibly modified) norm type \code{normtype}
            before last (internal) update}
      \item{risk}{(possibly [norm-]modified) risk \code{risk} from argument}
      \item{std}{(possibly modified) argument \code{std}}
      \item{iter}{number of iterations needed}
      \item{prec}{precision achieved}
      \item{b}{used clippng height \code{b}}
      \item{call}{call with which either \code{getLagrangeMultByIter} or
      \code{getLagrangeMultByOptim} was called }
}

\references{
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106-115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Ruckdeschel, P. and Rieder, H. (2004) Optimal Influence Curves for
  General Loss Functions. Statistics & Decisions \bold{22}: 201-223.

  Ruckdeschel, P. (2005) Optimally One-Sided Bounded Influence Curves.
  Mathematical Methods in Statistics \emph{14}(1), 105-131.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}.
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{InfRobModel-class}}}
%\examples{}
\concept{influence curve}
\keyword{robust}
