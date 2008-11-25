\name{outlyingPlotIC}
\alias{outlyingPlotIC}
\title{Function outlyingPlotIC in Package `RobAStBase' }
\description{outlyingPlotIC produces an outlyingness plot based on distances applied
to ICs}
\usage{
outlyingPlotIC(data, IC.x, IC.y = IC.x, dist.x = NormType(),
                         dist.y, cutoff.y = cutoff.chisq(), cutoff.x = cutoff.sememp(), ...,
                         cutoff.quantile.x = 0.95,
                         cutoff.quantile.y = cutoff.quantile.x,
                         id.n, lab.pts, adj, cex.idn,
                         col.idn, lty.cutoff, lwd.cutoff, col.cutoff,
                         main = gettext("Outlyingness by means of a distance-distance plot")
                         )
}
\arguments{
  \item{data}{data coercable to \code{matrix}; the data  at which to produce the \code{ddPlot}.}
  \item{IC.x}{object of class   \code{IC} the influence curve to produce
  the  distances for the \code{x} axis.}
  \item{IC.y}{object of class   \code{IC} the influence curve to produce
  the  distances for the \code{y} axis.}
  \item{\dots}{further arguments to be passed to \code{plot.default}, \code{text}, and \code{abline}}
  \item{dist.x}{object of class \code{NormType}; the distance for the \code{x} axis.}
  \item{dist.y}{object of class \code{NormType}; the distance for the \code{y} axis.}
  \item{cutoff.x}{object of class \code{cutoff}; the cutoff information for the \code{x} axis
  (the vertical line discriminating 'good' and 'bad' points).}
  \item{cutoff.y}{object of class \code{cutoff}; the cutoff information for the \code{y} axis
  (the horizontal line discriminating 'good' and 'bad' points).}
  \item{cutoff.quantile.x}{numeric; the cutoff quantile for the \code{x} axis.}
  \item{cutoff.quantile.y}{numeric; the cutoff quantile for the \code{y} axis.}
  \item{id.n}{a set of indices (or a corresponding logical vector); to select a subset
  of the data in argument \code{data}.}
  \item{lab.pts}{a vector of labels for the (unsubsetted) \code{data}.}
  \item{adj}{the corresponding argument for \code{\link[graphics]{text}} for
             labelling the outliers.}
  \item{cex.idn}{the corresponding \code{cex} argument for
                 \code{\link[graphics]{text}} for labelling the outliers.}
  \item{col.idn}{the corresponding \code{col} argument for
                 \code{\link[graphics]{text}} for labelling the outliers.}
  \item{lty.cutoff}{the corresponding \code{lty} argument for
                 \code{\link[graphics]{abline}} for drawing the cutoff lines.}
  \item{lwd.cutoff}{the corresponding \code{lwd} argument for
                 \code{\link[graphics]{abline}} for drawing the cutoff lines.}
  \item{col.cutoff}{the corresponding \code{col} argument for
                 \code{\link[graphics]{abline}} for drawing the cutoff lines.}
  \item{main}{the main title.}
}
\details{
calls a corresponding \code{\link{ddPlot}} method to produce the plot.
}
\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
  }

\value{
a list with items
\item{id.x}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{x}-cutoff}
\item{id.y}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{y}-cutoff}
\item{id.xy}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{x}-cutoff and the \code{y}-cutoff}
\item{qtx}{the quantiles of the distances of the (possibly transformed) data in \code{x} direction}
\item{qty}{the quantiles of the distances of the (possibly transformed) data in \code{y} direction}
\item{cutoff.x.v}{the cutoff value in \code{x} direction}
\item{cutoff.y.v}{the cutoff value in \code{y} direction}
}

\examples{
if(require(ROptEst)){
## generates normal location and scale family with mean = -2 and sd = 3
N0 <- NormLocationScaleFamily()
N0.IC0 <- optIC(model = N0, risk = asCov())
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))
N0.IC1 <- optIC(model = N0.Rob1, risk = asMSE())
xn <- c(rnorm(100),rcauchy(20)+20)
outlyingPlotIC(xn, IC.x=N0.IC0)
outlyingPlotIC(xn, IC.x=N0.IC1)
}
}
\keyword{hplot}

