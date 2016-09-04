\name{qqplot}
\docType{methods}
\title{Methods for Function qqplot in Package `RobAStBase'}
\usage{
qqplot(x, y, ...)
\S4method{qqplot}{ANY,RobModel}(x, y,
   n = length(x), withIdLine = TRUE, withConf = TRUE,
   withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ..., distance = NormType(),
    n.adj = TRUE)
\S4method{qqplot}{ANY,InfRobModel}(x, y, n = length(x), withIdLine = TRUE, 
withConf = TRUE, withConf.pw = withConf, withConf.sim = withConf,
  plot.it = TRUE, xlab = deparse(substitute(x)), ylab =
  deparse(substitute(y)), ..., cex.pts.fun = NULL, n.adj = TRUE)
\S4method{qqplot}{ANY,kStepEstimate}(x, y,
   n = length(x), withIdLine = TRUE, withConf = TRUE,
   withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...,
    exp.cex2.lbl = -.15,
    exp.cex2.pch = -.35,
    exp.fadcol.lbl = 1.85,
    exp.fadcol.pch = 1.85,
    bg = "white")
   }
\alias{qqplot}
\alias{qqplot-methods}
\alias{qqplot,ANY,RobModel-method}
\alias{qqplot,ANY,InfRobModel-method}
\alias{qqplot,ANY,kStepEstimate-method}

\arguments{
\item{x}{data to be checked for compatibility with distribution/model \code{y}.}
\item{y}{object of class \code{"RobModel"}, of class \code{"InfRobModel"} or of
class \code{"kStepEstimate"}.}
\item{n}{numeric; number of quantiles at which to do the comparison.}
\item{withIdLine}{logical; shall line \code{y = x} be plotted in?}
\item{withConf}{logical; shall confidence lines be plotted?}
\item{withConf.pw}{logical; shall pointwise confidence lines be plotted?}
\item{withConf.sim}{logical; shall simultaneous confidence lines be plotted?}
\item{plot.it}{logical; shall be plotted at all (inherited from \code{\link[stats:qqnorm]{qqplot}})?}
\item{xlab}{x-label}
\item{ylab}{y-label}
\item{\dots}{further parameters for method \code{qqplot} with signature
\code{ANY,ProbFamily} (see \code{\link[distrMod]{qqplot}}) or with function \code{plot}}
\item{cex.pts.fun}{rescaling function for the size of the points to be plotted;
        either \code{NULL} (default), then \code{log(1+abs(x))} is used,
        or a function which is then used. }
\item{n.adj}{logical; shall sample size be adjusted for possible outliers according
to radius of the corresponding neighborhood?}
\item{distance}{a function mapping observations \code{x} to the positive reals;
used to determine the size of the plotted points (the larger \code{distance(x)},
the smaller the points are plotted.}
\item{exp.cex2.lbl}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to magnify the labels.}
\item{exp.cex2.pch}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to magnify the symbols.}
\item{exp.fadcol.lbl}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to find out-fading colors.}
\item{exp.fadcol.pch}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to find out-fading colors.}
\item{bg}{background color to fade against}
}
\description{
  We generalize function \code{\link[stats:qqnorm]{qqplot}} from package \pkg{stats} to
  be applicable to distribution and probability model objects. In this context,
  \code{qqplot} produces a QQ plot of data (argument \code{x}) against
   a (model) distribution. For arguments \code{y} of class \code{RobModel},
   points at a high \dQuote{distance} to the model
   are plotted smaller. For arguments \code{y} of class \code{kStepEstimate},
   points at with low weight in the [p]IC are plotted bigger and their
   color gets faded out slowly.
   Graphical parameters may be given as arguments to \code{qqplot}.
}
\value{
    As for function \code{\link[stats:qqnorm]{qqplot}} from package \pkg{stats}: a
  list with components
\item{x}{The x coordinates of the points that were/would be plotted}
  \item{y}{The corresponding quantiles of the second distribution,
           \emph{including \code{\link{NA}}s}.}
}
\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language}.
  Wadsworth & Brooks/Cole.
}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
}
\seealso{
  \code{\link[stats:qqnorm]{qqplot}} from package \pkg{stats} -- the standard QQ plot
  function,  \code{\link[distr]{qqplot}} from package \pkg{distr} for
  comparisons of distributions, and
 \code{\link[distrMod]{qqplot}} from package \pkg{distrMod} (which
  is called intermediately by this method), as well as
  \code{\link{qqbounds}}, used by \code{qqplot} to produce confidence
   intervals.
}
\details{
\describe{
\item{qqplot}{\code{signature(x = "ANY", y = "RobModel")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of distribution of robust model \code{y}.}
\item{qqplot}{\code{signature(x = "ANY", y = "InfRobModel")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of distribution of infinitesimally robust model \code{y}.}
\item{qqplot}{\code{signature(x = "ANY", y = "kStepEstimate")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of the model distribution of model at which
the corresponding \code{kStepEstimate} \code{y} had been calibrated at.
By default, if the [p]IC of the \code{kStepEstimate} is of class
\code{HampIC}, i.e.; has a corresponding weight function,
points (and, if \code{withLab==TRUE}, labels) are
scaled and faded according to this weight function. Corresponding
arguments \code{exp.cex2.pch} and \code{exp.fadcol.pch} control this
scaling and fading, respectively
(and analogously \code{exp.cex2.lbl} and \code{exp.fadcol.lbl} for the labels).
The choice of these arguments has to be done on a case-by-case basis.
Positive exponents induce fading, magnification with increasing weight,
for negative exponents the same is true for decreasing weight; higher
(absolute) values increase the speed of fading / magnification.
}
}
}

\examples{
qqplot(rnorm(40, mean = 15, sd = sqrt(30)), Chisq(df=15))
RobM <- InfRobModel(center = NormLocationFamily(mean=13,sd=sqrt(28)),
                    neighbor = ContNeighborhood(radius = 0.4))
x <- rnorm(20, mean = 15, sd = sqrt(30))
qqplot(x, RobM)
\donttest{ 
## \donttest to reduce check time
qqplot(x, RobM, alpha.CI=0.9, add.points.CI=FALSE)
}
## further examples for ANY,kStepEstimator-method
## in example to roptest() in package ROptEst
}
\keyword{hplot}
\keyword{distribution}
