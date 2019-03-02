\name{returnlevelplot}
\docType{methods}
\title{Methods for Function returnlevelplot in Package `RobAStBase'}
\usage{
returnlevelplot(x, y, ...)
\S4method{returnlevelplot}{ANY,RobModel}(x, y,
   n = length(x), withIdLine = TRUE, withConf = TRUE,
   withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ..., distance = NormType(),
    n.adj = TRUE)
\S4method{returnlevelplot}{ANY,InfRobModel}(x, y, n = length(x), withIdLine = TRUE,
withConf = TRUE, withConf.pw = withConf, withConf.sim = withConf,
  plot.it = TRUE, xlab = deparse(substitute(x)), ylab =
  deparse(substitute(y)), ..., cex.pts.fun = NULL, n.adj = TRUE)
\S4method{returnlevelplot}{ANY,kStepEstimate}(x, y,
   n = length(x), withIdLine = TRUE, withConf = TRUE,
   withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...,
    exp.cex2.lbs = -.15,
    exp.cex2.pts = -.35,
    exp.fadcol.lbs = 1.85,
    exp.fadcol.pts = 1.85,
    bg = "white")
   }
\alias{returnlevelplot}
\alias{returnlevelplot-methods}
\alias{returnlevelplot,ANY,RobModel-method}
\alias{returnlevelplot,ANY,InfRobModel-method}
\alias{returnlevelplot,ANY,kStepEstimate-method}

\arguments{
\item{x}{data to be checked for compatibility with distribution/model \code{y}.}
\item{y}{object of class \code{"RobModel"}, of class \code{"InfRobModel"} or of
class \code{"kStepEstimate"}.}
\item{n}{numeric; number of quantiles at which to do the comparison.}
\item{withIdLine}{logical; shall line \code{y = x} be plotted in?}
\item{withConf}{logical; shall confidence lines be plotted?}
\item{withConf.pw}{logical; shall pointwise confidence lines be plotted?}
\item{withConf.sim}{logical; shall simultaneous confidence lines be plotted?}
\item{plot.it}{logical; shall be plotted at all (inherited from
\code{\link[distrMod:returnlevelplot]{returnlevelplot}})?}
\item{xlab}{x-label}
\item{ylab}{y-label}
\item{\dots}{further parameters for method \code{returnlevelplot} with signature
\code{ANY,ProbFamily} (see \code{\link[distrMod]{returnlevelplot}}) or with function
 \code{plot}}
\item{cex.pts.fun}{rescaling function for the size of the points to be plotted;
        either \code{NULL} (default), then \code{log(1+abs(x))} is used,
        or a function which is then used. }
\item{n.adj}{logical; shall sample size be adjusted for possible outliers according
to radius of the corresponding neighborhood?}
\item{distance}{a function mapping observations \code{x} to the positive reals;
used to determine the size of the plotted points (the larger \code{distance(x)},
the smaller the points are plotted.}
\item{exp.cex2.lbs}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to magnify the labels.}
\item{exp.cex2.pts}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to magnify the symbols.}
\item{exp.fadcol.lbs}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to find out-fading colors.}
\item{exp.fadcol.pts}{for objects \code{kStepEstimate} based on a [p]IC of class \code{HampIC}:
exponent for the weights of this [p]IC used to find out-fading colors.}
\item{bg}{background color to fade against}
}
\description{
  We generalize function \code{\link[distrMod]{returnlevelplot}} from package \pkg{distrMod} to
  be applicable to distribution and probability model objects. In this context,
  \code{returnlevelplot} produces a rescaled QQ plot of data (argument \code{x})
  against a (model) distribution. For arguments \code{y} of class \code{RobModel},
   points at a high \dQuote{distance} to the model
   are plotted smaller. For arguments \code{y} of class \code{kStepEstimate},
   points at with low weight in the [p]IC are plotted bigger and their
   color gets faded out slowly. This parallels the behaviour of the respective
   \code{qqplot} methods.
   Graphical parameters may be given as arguments to \code{returnlevelplot}.
}
\value{
    As for function \code{\link[distrMod:returnlevelplot]{returnlevelplot}} from package \pkg{stats}.
}
\note{
The confidence bands given in our version of the return level plot differ
 from the ones given in package \pkg{ismev}. We use non-parametric bands,
 hence also allow for non-parametric deviances from the model, whereas in
 in package \pkg{ismev} they are based on profiling, hence only check for
 variability within the parametric class.
}
\references{
  ismev: An   Introduction to Statistical Modeling of Extreme Values. R package
  version 1.39. https://CRAN.R-project.org/package=ismev; original S functions
  written by Janet E. Heffernan with R port and R documentation provided by
  Alec G. Stephenson. (2012).

  Coles, S. (2001). \emph{An introduction
  to statistical modeling of extreme values.} London: Springer.
}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
}
\seealso{
  \code{\link[stats:qqnorm]{qqplot}} from package \pkg{stats} -- the standard QQ plot
  function,  \code{\link[distrMod]{returnlevelplot}} from package \pkg{distrMod} (which
  is called intermediately by this method), as well as
  \code{\link{qqbounds}}, used by \code{returnlevelplot} to produce confidence
   intervals.
}
\details{
\describe{
\item{returnlevelplot}{\code{signature(x = "ANY", y = "RobModel")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of distribution of robust model \code{y}.}
\item{returnlevelplot}{\code{signature(x = "ANY", y = "InfRobModel")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of distribution of infinitesimally robust model \code{y}.}
\item{returnlevelplot}{\code{signature(x = "ANY", y = "kStepEstimate")}:
produces a QQ plot of a dataset \code{x} against the theoretical
quantiles of the model distribution of model at which
the corresponding \code{kStepEstimate} \code{y} had been calibrated at.
By default, if the [p]IC of the \code{kStepEstimate} is of class
\code{HampIC}, i.e.; has a corresponding weight function,
points (and, if \code{withLab==TRUE}, labels) are
scaled and faded according to this weight function. Corresponding
arguments \code{exp.cex2.pts} and \code{exp.fadcol.pts} control this
scaling and fading, respectively
(and analogously \code{exp.cex2.lbs} and \code{exp.fadcol.lbs} for the labels).
The choice of these arguments has to be done on a case-by-case basis.
Positive exponents induce fading, magnification with increasing weight,
for negative exponents the same is true for decreasing weight; higher
(absolute) values increase the speed of fading / magnification.
}
}
}

\examples{
returnlevelplot(rnorm(40, mean = 15, sd = sqrt(30)), Chisq(df=15))
RobM <- InfRobModel(center = NormLocationFamily(mean=13,sd=sqrt(28)),
                    neighbor = ContNeighborhood(radius = 0.4))
\donttest{
## \donttest to reduce check time
x <- rnorm(20, mean = 15, sd = sqrt(30))
returnlevelplot(x, RobM)
returnlevelplot(x, RobM, alpha.CI=0.9, add.points.CI=FALSE)
}
## further examples for ANY,kStepEstimator-method
## in example to roptest() in package ROptEst
}
\keyword{hplot}
\keyword{distribution}
