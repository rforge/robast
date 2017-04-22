\name{QuantileBCCEstimator}
\alias{QuantileBCCEstimator}
\alias{.QBCC}

\title{ Function to compute QuantileBCC estimates for the Weibull Family}
\description{
  Function \code{QuantileBCCEstimator} computes QuantileBCC estimator
  (for the Weibull) at real data and returns an object of class \code{Estimate}.
}
\usage{
QuantileBCCEstimator(x, p1 = 1/3, p2 = 2/3,
            name, Infos, nuis.idx = NULL,
            trafo = NULL, fixed = NULL, na.rm = TRUE,
            ...)
.QBCC(x, p1 = 1/3, p2 = 2/3)
}
\arguments{
  \item{x}{ (empirical) data }
  \item{p1,p2}{levels of the quantiles; maximal breakdown point is achieved
               for \eqn{p1=p2-p1=1-p2=1/3} which is the default.}
  \item{name}{ optional name for estimator. }
  \item{Infos}{ character: optional informations about estimator }
  \item{nuis.idx}{ optionally the indices of the estimate belonging
                  to nuisance parameter}
  \item{fixed}{ optionally (numeric) the fixed part of the parameter}
  \item{trafo}{ an object of class \code{MatrixorFunction} -- a transformation
  for the main parameter}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{\dots}{not yet used. }
}
\details{
  The actual work is done in \code{.QBCC}.
  The wrapper \code{QuantileBCCEstimator} pre-treats the data,
  and constructs a respective \code{Estimate} object.
}
\value{
  \item{.QuantileBCCEstimator}{A numeric vector of length \code{2} with components
   named \code{scale} and \code{shape}. }
  \item{QuantileBCCEstimator}{An object of S4-class \code{"Estimate"}. }
}
\references{
Boudt, K., Caliskan, D., Croux, C. (2011): Robust explicit estimators of Weibull parameters.
\emph{Metrika}, \bold{73} (2), 187-209.
}

%\references{  }
\author{Nataliya Horbenko \email{nhorbenko@gmail.com},\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{ParamFamily-class}}, \code{\link{ParamFamily}}, 
         \code{\link{Estimate-class}} }
\examples{
## (empirical) Data
set.seed(123)
distroptions("withgaps"=FALSE)
x <- rweibull(50, scale = 0.5, shape = 3)
##
QuantileBCCEstimator(x = x)
}
\keyword{univar}
