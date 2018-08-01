\name{optIC}
\alias{optIC}
\alias{optIC-methods}
\alias{optIC,InfRobModel,asRisk-method}
\alias{optIC,InfRobModel,asUnOvShoot-method}
\alias{optIC,FixRobModel,fiUnOvShoot-method}

\title{Generic function for the computation of optimally robust ICs}
\description{
  Generic function for the computation of optimally robust ICs.
}
\usage{
optIC(model, risk, ...)

\S4method{optIC}{InfRobModel,asRisk}(model, risk, z.start = NULL, A.start = NULL,
                                     upper = 1e4, lower = 1e-4,
                                     OptOrIter = "iterate", maxiter = 50,
                                     tol = .Machine$double.eps^0.4, warn = TRUE, 
                                     noLow = FALSE, verbose = NULL, ...,
                                     .withEvalAsVar = TRUE, withMakeIC = FALSE,
                                     returnNAifProblem = FALSE, modifyICwarn = NULL)

\S4method{optIC}{InfRobModel,asUnOvShoot}(model, risk, upper = 1e4,
                                          lower = 1e-4, maxiter = 50,
                                          tol = .Machine$double.eps^0.4,
                                          withMakeIC = FALSE, warn = TRUE,
                                          verbose = NULL, modifyICwarn = NULL)

\S4method{optIC}{FixRobModel,fiUnOvShoot}(model, risk, sampleSize, upper = 1e4, lower = 1e-4,
                                          maxiter = 50, tol = .Machine$double.eps^0.4, 
                                          withMakeIC = FALSE, warn = TRUE,
                                          Algo = "A", cont = "left",
                                          verbose = NULL, modifyICwarn = NULL)
}
\arguments{
  \item{model}{ probability model. }
  \item{risk}{ object of class \code{"RiskType"}. }
  \item{\dots}{ additional parameters. }
  \item{z.start}{ initial value for the centering constant. }
  \item{A.start}{ initial value for the standardizing matrix. }
  \item{upper}{ upper bound for the optimal clipping bound. }
  \item{lower}{ lower bound for the optimal clipping bound. }
  \item{maxiter}{ the maximum number of iterations. }
  \item{tol}{ the desired accuracy (convergence tolerance).}
  \item{warn}{ logical: print warnings. }  
  \item{sampleSize}{ integer: sample size. }
  \item{Algo}{ "A" or "B". }
  \item{cont}{ "left" or "right". }
  \item{noLow}{ logical: is lower case to be computed? }
  \item{OptOrIter}{character; which method to be used for determining Lagrange
  multipliers \code{A} and \code{a}: if (partially) matched to \code{"optimize"},
  \code{getLagrangeMultByOptim} is used; otherwise: by default, or if matched to
  \code{"iterate"} or to \code{"doubleiterate"},
  \code{getLagrangeMultByIter} is used. More specifically,
  when using \code{getLagrangeMultByIter}, and if argument \code{risk} is of
  class \code{"asGRisk"}, by default and if matched to \code{"iterate"}
  we use only one (inner) iteration, if matched to \code{"doubleiterate"}
  we use up to \code{Maxiter} (inner) iterations. }
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed. }
  \item{.withEvalAsVar}{logical (of length 1):
     if \code{TRUE}, risks based on covariances are to be
     evaluated (default), otherwise just a call is returned. }
  \item{withMakeIC}{logical; if \code{TRUE} the [p]IC is passed through
        \code{makeIC} before return.}
  \item{returnNAifProblem}{logical (of length 1):
     if \code{TRUE} (not the default), in case of convergence problems in
     the algorithm, returns \code{NA}. }
  \item{modifyICwarn}{logical: should a (warning) information be added if
     \code{modifyIC} is applied and hence some optimality information could
     no longer be valid? Defaults to \code{NULL} in which case this value
     is taken from \code{RobAStBaseOptions}.}
}
\details{ In case of the finite-sample risk \code{"fiUnOvShoot"} one can choose
  between two algorithms for the computation of this risk where the least favorable
  contamination is assumed to be left or right of some bound. For more details
  we refer to Section 11.3 of Kohl (2005).}
\value{Some optimally robust IC is computed.}
\section{Methods}{
\describe{
  \item{model = "InfRobModel", risk = "asRisk"}{ 
    computes optimally robust influence curve for 
    robust models with infinitesimal neighborhoods and
    various asymptotic risks. }

  \item{model = "InfRobModel", risk = "asUnOvShoot"}{ 
    computes optimally robust influence curve for 
    robust models with infinitesimal neighborhoods and
    asymptotic under-/overshoot risk. }

  \item{model = "FixRobModel", risk = "fiUnOvShoot"}{ 
    computes optimally robust influence curve for 
    robust models with fixed neighborhoods and
    finite-sample under-/overshoot risk. }
}}
\references{
 
  Huber, P.J. (1968) Robust Confidence Limits. Z. Wahrscheinlichkeitstheor.
  Verw. Geb. \bold{10}:269--278.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Kohl, M. and Ruckdeschel, P. (2010): R package distrMod: 
  Object-Oriented Implementation of Probability Models. 
  J. Statist. Softw.  \bold{35}(10), 1--27
  
  Kohl, M. and Ruckdeschel, P., and Rieder, H. (2010):
  Infinitesimally Robust Estimation in General Smoothly Parametrized Models.
  \emph{Stat. Methods Appl.}, \bold{19}, 333--354. 
  
  Rieder, H. (1980) Estimates derived from robust tests. Ann. Stats. \bold{8}: 106--115.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \bold{17}(1) 13-40.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[RobAStBase]{InfluenceCurve-class}}, \code{\link[distrMod]{RiskType-class}}}
\examples{
B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)
}
\concept{robust influence curve}
\concept{influence curve}
\keyword{robust}
