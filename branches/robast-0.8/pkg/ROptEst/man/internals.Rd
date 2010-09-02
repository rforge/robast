\name{internals_for_ROptEst}
\alias{internals_for_ROptEst}
\alias{.checkUpLow}
\alias{.getUpperSol}
\alias{.getLowerSol}
\alias{.getLowUpB}
\alias{.checkPIC}
\alias{.LowerCaseMultivariate}
\alias{.LowerCaseMultivariateTV}
\alias{.getSB}

\title{Internal / Helper functions of package ROptEst}

\description{
These functions are used internally by package \pkg{ROptEst}.}

\usage{
### helper function to check whether given b is in (bmin, bmax)
###        if not returns corresponding upper / lower case solution

.checkUpLow(L2deriv, b, risk, neighbor, biastype, normtype,
                        Distr, Finfo, DistrSymm, L2derivSymm,
                        L2derivDistrSymm, z.start, A.start, trafo, maxiter,
                        tol, QuadForm, verbose, nrvalpts, warn)
                        
### helper function to return the upper case solution if r=0
.getUpperSol(L2deriv, b, radius, risk, neighbor, biastype,
                       normtype, Distr, Finfo, trafo,
                       QuadForm, verbose, warn)

### helper function to return the lower case solution if b-search was not successful
.getLowerSol(L2deriv, risk, neighbor, Distr, DistrSymm,
                         L2derivSymm, L2derivDistrSymm,
                         z.start, A.start, trafo,
                         maxiter, tol, warn, Finfo, QuadForm, verbose)


### helper function to return upper & lower bounds for b for b-search
.getLowUpB(L2deriv, Finfo, Distr, normtype, z, A, radius, iter)

### helper function to check whether (TotalVariation) weight w has already been modified
.isVirginW(w)

### helper function to check whether (intermediate) results give a pIC
.checkPIC(L2deriv, neighbor, Distr, trafo, z, A, w, z.comp, A.comp, ...)

.LowerCaseMultivariate(L2deriv, neighbor, biastype,
             normtype, Distr, Finfo, trafo, z.start = NULL,
             A.start = NULL, z.comp = NULL, A.comp = NULL,
             maxiter, tol, verbose = NULL)

.LowerCaseMultivariateTV(L2deriv, neighbor, biastype,
             normtype, Distr, Finfo, trafo,
             A.start,  maxiter, tol,
             verbose = NULL)

.getSB(IC,neighbor)
}

\arguments{
  \item{L2deriv}{ L2-derivative of some L2-differentiable family
    of probability measures. }
  \item{b}{numeric; clipping bound under consideration.}
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{neighbor}{ object of class \code{"Neighborhood"}. }
  \item{biastype}{object of class \code{"BiasType"} --- the bias type with we work.}
  \item{normtype}{object of class \code{"NormType"} --- the norm type with we work.}
  \item{Distr}{ object of class \code{"Distribution"}. }
  \item{Finfo}{ Fisher information matrix. }
  \item{DistrSymm}{ object of class \code{"DistributionSymmetry"}. }
  \item{L2derivSymm}{ object of class \code{"FunSymmList"}. }
  \item{L2derivDistrSymm}{ object of class \code{"DistrSymmList"}. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{trafo}{ matrix: transformation of the parameter. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{QuadForm}{ matrix of (or which may coerced to) class
                  \code{PosSemDefSymmMatrix} for use of different
                  (standardizing) norm }
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed. }
  \item{nrvalpts}{integer: number of evaluation points.}
  \item{warn}{ logical: print warnings. }
  \item{radius}{radius of the neighborhood.}
  \item{z}{centering constant (in \code{k}-space)}
  \item{A}{standardizing matrix.}
  \item{w}{a weight of class \code{"BdStWeight"}}
  \item{z.comp}{logical vector: indicator which components of \code{z} need
                to be computed}
  \item{A.comp}{logical matrix: indicator which components of \code{A} need
                to be computed}
  \item{iter}{the number of iterations computed so far; used for specifying
              a different value of the clipping component of the weight in
              total variation case in the very first iteration.}
  \item{IC}{some IC of class \code{"IC"}}
  \item{\dots}{further arguments to be passed on \code{E()}.}
}

\details{
\code{.checkUpLow} checks whether the given clipping height \code{b} lies in
    \eqn{(b_{\rm\scriptstyle min},b_{\rm\scriptstyle min})}{(b_min,b_max)};
\code{.getUpperSol} determines the upper case/classical solution and computes
  corresponding risks
\code{.getLowerSol} determines the lower case (minimax bias) solution and computes
  corresponding risks
\code{.getLowUpB} determines a search interval for \code{b} to given radius
\code{r}, i.e., lower and upper bounds for
\eqn{(b_{\rm\scriptstyle min},b_{\rm\scriptstyle min})}{(b_min,b_max)}
\code{.isVirginW} checks whether the (total variation) weight \code{w} in
the argument has already been modified since creation (\code{TRUE} if not)
\code{.checkPIC} checks whether (intermediate) results give a pIC
\code{.LowerCaseMultivariatefunction} determines the Lagrange multipliers for
the multivariate lower case solution for convex contamination
by solving a corresponding dual problem (Rieder[94],p.199 eq.(18)).
\code{.LowerCaseMultivariatefunctionTV} determines the Lagrange multipliers for
the multivariate lower case solution for total variation in dimension \eqn{p=1}
and \eqn{k>1} by solving a corresponding dual problem (Rieder[94],p.205 eq.(58)).
\code{.getSB} computes the bias and (the square root of the trace of) the variance
 of the IC.
}


\value{
\item{.checkUpLow}{a list with items \code{up} (logical)
                  \code{low} (logical) and \code{res} a return list for
                  \code{getInfRobIC} or \code{NULL}: if \code{up} is \code{TRUE},
                  \code{res} contains the upper case solution,
                  if \code{low} is \code{TRUE} the lower case solution
                  and if both \code{low} and \code{up} are \code{FALSE}
                  \code{NULL}.}
\item{.getUpperSol}{a return list for \code{getInfRobIC}}
\item{.getLowerSol}{a return list for \code{getInfRobIC}}
\item{.checkUpLow}{a list with items \code{lower} and \code{upper} (both numeric).}
\item{.isVirginW}{\code{TRUE} or \code{FALSE}}
\item{.checkPIC}{nothing is returned; precision values are issued.}
\item{.LowerCaseMultivariatefunction}{a list with elements
\code{erg} the return value of the corresponding call to \code{optim},
\code{w} the corresponding weight (of class \code{HampelWeight}),
\code{normtype} the norm ---possibly modified,
\code{z.comp} the indices of non-null entries in \code{z}.}
\item{.LowerCaseMultivariatefunctionTV}{a list with elements
\code{A} the standardization matrix,
\code{b} the minimax bias,
\code{w} the corresponding weight (of class \code{BdStWeight}),
\code{a} the corresponding infimum of the [p]IC.}
\item{.getSB}{a list with elements \code{s}, (the square root of the trace of) 
the variance and \code{b}, the bias.}
}

\references{
  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
}

\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
  }

\keyword{internal}
\concept{utilities}
