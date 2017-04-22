\name{PickandsEstimator}
\alias{PickandsEstimator}
\alias{.PickandsEstimator}

\title{ Function to compute Pickands estimates for the GPD and GEVD}
\description{
  Function \code{PickandsEstimator} computes Pickands estimator
  (for the GPD and GEVD) at real data and returns an object of class \code{Estimate}.
}
\usage{
PickandsEstimator(x, ParamFamily=GParetoFamily(), alpha=2,
            name, Infos, nuis.idx = NULL,
            trafo = NULL, fixed = NULL, na.rm = TRUE,
            ...)
.PickandsEstimator(x, alpha=2, GPD.l = TRUE)
}
\arguments{
  \item{x}{ (empirical) data }
  \item{alpha}{ numeric \eqn{> 1}; determines the variant of the Pickands-Estimator
   based on matching the empirical quantiles  to levels
   \eqn{a_1=1-1/\alpha}{\code{a1=1-1/alpha}} and
   \eqn{a_2=1-1/\alpha^2}{\code{a2=1-1/alpha^2}} (in the GPD case) resp.
   \eqn{a_1=\exp(-1/\alpha)}{\code{a1=exp(-1/alpha)}} and
   \eqn{a_1=\exp(-1/\alpha^2)}{\code{a2=exp(-1/alpha^2)}} (in the GEVD case)
   against the population counter parts. The ''classical'' Pickands Estimator
   building up on the median is obtained for \code{alpha=2} for the GPD
   and for \code{alpha = 1/log(2)} for the GEVD.
   If \code{alpha} is missing we set it to the optimal value (see note below).}
  \item{ParamFamily}{an object of class \code{"GParetoFamily"} or
                     \code{"GEVFamily"}. }
  \item{name}{ optional name for estimator. }
  \item{Infos}{ character: optional informations about estimator }
  \item{nuis.idx}{ optionally the indices of the estimate belonging
                  to nuisance parameter}
  \item{fixed}{ optionally (numeric) the fixed part of the parameter}
  \item{trafo}{ an object of class \code{MatrixorFunction} -- a transformation
  for the main parameter}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{\dots}{not yet used. }
  \item{GPD.l}{logical: if  \code{TRUE} the variant for GPD is used, else for
   GEVD.}
}
\details{
  The actual work is done in \code{.PickandsEstimator}.
  The wrapper \code{PickandsEstimator} pre-treats the data,
  and constructs a respective \code{Estimate} object.
}
\note{
The scale estimate we use, i.e., with scale = \eqn{\beta}{beta}
and shape = \eqn{\xi}{xi}, we estimate scale by
\eqn{\beta= \xi a_1/(\alpha^\xi-1)}{\code{beta= xi*a1/(alpha^xi-1)}},  differs from
the one given in the original reference, where it was
\eqn{\beta= \xi a_1^2/(a_2-2a_1)}{\code{beta= xi * a1^2 /(a2-2*a1)}}.
The one chosen here avoids taking differences \eqn{a_2-2a_1}{a2-2*a1}
hence does not require \eqn{a_2 > 2a_1}{a2>2*a1}; this leads to
(functional) breakdown point (bdp)
\deqn{\min(a_1,1-a_2,a_2-a_1)}{min(a1,1-a2,a2-a1)}
which is independent \eqn{\xi}{xi}, whereas the original setting leads to
a bdp  which is depending on \eqn{\xi}{xi}
\deqn{\min(a_1,1-a_2,a_2-1+(2\alpha^\xi-1)^{-1/\xi})\qquad
\mbox{for GPD}}{\code{min(a1,1-a2,a2-1+(2*alpha^xi-1)^(-1/xi))}   for GPD}
\deqn{\min(a_1,1-a_2,a_2-\exp(-(2\alpha^\xi-1)^{-1/\xi})) \qquad
\mbox{for GEVD}}{\code{min(a1,1-a2,a2-exp(-(2*alpha^xi-1)^(-1/xi)))}
for GEVD}. As a consequence our setting, the bdp-optimal choice of
\eqn{\alpha}{alpha} for GDP is \eqn{2} leading to bdp \eqn{1/4}, and
\eqn{2.248} for GEVD leading to bdp \eqn{0.180}. For comparison, with the
original setting, at \eqn{\xi=0.7}{xi=0.7}, this  gives optimal bdp's
\eqn{0.070} and \eqn{0.060} for GPD and GEVD, respectively.
The standard choice of \eqn{\alpha}{alpha} such that \eqn{a_1}{a1}
gives the median (\eqn{\alpha=2}{alpha=2} in the GPD and
\eqn{\alpha=1/\log(2)}{alpha=1/log(2)} in the GEVD) in our setting gives
bdp's of \eqn{1/4} and \eqn{0.119} for GPD and GEVD, respectively, and
in the original setting, at \eqn{\xi=0.7}{xi=0.7}, gives bdp's
\eqn{0.064} and \eqn{0.023}.
}
\value{
  \item{.PickandsEstimator}{A numeric vector of length \code{2} with components
   named \code{scale} and \code{shape}. }
  \item{PickandsEstimator}{An object of S4-class \code{"Estimate"}. }
}
\references{
P. Ruckdeschel, N. Horbenko (2012): Yet another breakdown point notion:
EFSBP --illustrated at scale-shape models. \emph{Metrika}, \bold{75}(8),
1025--1047.\cr

J. Pickands (1975): Statistical inference using extreme order statistics.
\emph{Ann. Stat.} \bold{3}(1), 119--131.
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
x <- rgpd(50, scale = 0.5, shape = 3)
y <- rgev(50, scale = 0.5, shape = 3)
## parametric family of probability measures
P <- GParetoFamily(scale = 1, shape = 2)
G <- GEVFamily(scale = 1, shape = 2)
##
PickandsEstimator(x = x, ParamFamily = P)
PickandsEstimator(x = y, ParamFamily = G)
}
\keyword{univar}
