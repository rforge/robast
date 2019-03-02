\name{GEVFamily}
\alias{GEVFamily}

\title{Generating function for families of Generalized Extreme Value distributions }
\description{
  Generates an object of class \code{"GEVFamily"} which
  represents a Generalized EV family.
}
\usage{
GEVFamily(loc = 0, scale = 1, shape = 0.5, of.interest = c("scale", "shape"),
          p = NULL, N = NULL, trafo = NULL, start0Est = NULL, withPos = TRUE,
          secLevel = 0.7, withCentL2 = FALSE, withL2derivDistr  = FALSE,
          withMDE = FALSE, ..ignoreTrafo = FALSE, ..withWarningGEV = TRUE)
}
\arguments{
  \item{loc}{ real: known/fixed threshold/location parameter }
  \item{scale}{ positive real: scale parameter }
  \item{shape}{ positive real: shape parameter }
  \item{of.interest}{ character: which parameters, transformations are of interest.\cr
              possibilites are: "scale", "shape", "quantile", "expected loss",
              "expected shortfall"; a maximum number of two of these may be selected }
  \item{p}{real or NULL: probability needed for quantile and expected shortfall }
  \item{N}{real or NULL: expected frequency for expected loss }
  \item{trafo}{ matrix or NULL: transformation of the parameter }
  \item{start0Est}{ startEstimator --- if \code{NULL} \code{\link{PickandsEstimator}} is used }
  \item{withPos}{ logical of length 1: Is shape restricted to positive values? }
  \item{secLevel}{ a numeric of length 1:
       In the ideal GEV model, for each observastion \eqn{X_i}{Xi}, the expression
       \eqn{1+\frac{{\rm shape}(X_i-{\rm loc})}{{\rm scale}}}{1+shape(Xi-loc)/scale}
       must be positive, which in principle could be attacked by a single outlier.
       Hence for sample size \eqn{n} we allow for \eqn{\varepsilon n}{eps n}
       violations, interpreting the violations as outliers. Here
       \eqn{\varepsilon = {\tt secLevel}/\sqrt{n}}{eps = secLevel/sqrt(n)}. }
  \item{withCentL2}{logical: shall L2 derivative be centered by substracting
       the E()? Defaults to \code{FALSE}, but higher accuracy can be achieved
       when set to \code{TRUE}.}
  \item{withL2derivDistr}{logical: shall the distribution of the L2 derivative
      be computed? Defaults to \code{FALSE} (to speed up computations).}
  \item{withMDE}{logical: should Minimum Distance Estimators be used to
                 find a good starting value for the parameter search?
                 Defaults to \code{FALSE}  (to speed up computations).
                 We have seen cases though, where the use of the then
                 employed \code{PickandsEstimator} was drastically misleading
                 and subsequently led to bad estimates where it is used
                 as starting value; so where feasible it is a good idea
                 to also try argument \code{withMDE=TRUE} for control purposes.}
  \item{..ignoreTrafo}{logical: only used internally in \code{kStepEstimator}; do not change this.}
  \item{..withWarningGEV}{logical: shall warnings be issued if shape is large?}
}
\details{
  The slots of the corresponding L2 differentiable 
  parameteric family are filled.
}
\value{Object of class \code{"GEVFamily"}}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to 
  the Asymptotic Theory of Robustness}. Bayreuth: Dissertation.\cr

M.~Kohl, P. Ruckdeschel, H.~Rieder (2010):
Infinitesimally Robust Estimation in General Smoothly Parametrized Models.
\emph{Stat. Methods Appl.}, \bold{19}, 333--354.\cr

Ruckdeschel, P. and Horbenko, N. (2012): Yet another breakdown point notion:
EFSBP --illustrated at scale-shape models. \emph{Metrika}, \bold{75}(8),
1025--1047.
  }

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}\cr
        Nataliya Horbenko \email{nhorbenko@gmail.com}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\linkS4class{GPareto}}}
\examples{
(G1 <- GEVFamily())
FisherInfo(G1)
checkL2deriv(G1)
}
\concept{Generalized Pareto model}
\keyword{models}
