\name{outlyingPlotIC}
\alias{outlyingPlotIC}
\title{Function outlyingPlotIC in Package `RobAStBase' }
\description{outlyingPlotIC produces an outlyingness plot based on distances applied
to ICs}
\usage{
outlyingPlotIC(data,IC.x, IC.y = IC.x, dist.x = NormType(), dist.y, 
 cutoff.x = cutoff.sememp(0.95), cutoff.y = cutoff.chisq(0.95),  ...,
 cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
 id.n, cex.pts = 1, lab.pts, jitter.pts = 0, alpha.trsp = NA, adj, cex.idn,
 col.idn,  lty.cutoff,  lwd.cutoff,  col.cutoff, text.abline = TRUE,
 text.abline.x = NULL, text.abline.y = NULL, cex.abline = par("cex"),
 col.abline = col.cutoff, font.abline = par("font"), adj.abline = c(0,0),
 text.abline.x.x = NULL, text.abline.x.y = NULL, text.abline.y.x = NULL,
 text.abline.y.y = NULL, text.abline.x.fmt.cx = "\%7.2f",
 text.abline.x.fmt.qx = "\%4.2f\%\%", text.abline.y.fmt.cy = "\%7.2f",
 text.abline.y.fmt.qy = "\%4.2f\%\%", robCov.x = TRUE, robCov.y = TRUE,
 tf.x = NULL,tf.y = NULL, jitter.fac=10, jitter.tol=.Machine$double.eps,
 doplot = TRUE,
 main = gettext("Outlyingness \n by means of a distance-distance plot")
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
  \item{cex.pts}{the corresponding \code{cex} argument for plotted points.}
  \item{lab.pts}{a vector of labels for the (unsubsetted) \code{data}.}
  \item{jitter.pts}{the corresponding \code{jitter} argument for plotted points;
    may be a vector of length 2 -- for separate factors for x- and y-coordinate.}
  \item{alpha.trsp}{alpha transparency to be added ex post to colors
        \code{col.pch} and \code{col.lbl}; if one-dim and NA all colors are
        left unchanged. Otherwise, with usual recycling rules \code{alpha.trsp}
        gets shorted/prolongated to length the data-symbols to be plotted.
        Coordinates of this vector \code{alpha.trsp} with NA are left unchanged,
        while for the remaining ones, the alpha channel in rgb space is set
        to the respective coordinate value of \code{alpha.trsp}. The non-NA
        entries must be integers in [0,255] (0 invisible, 255 opaque).}
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
  \item{text.abline}{vector of logicals (cast to length 2): shall text be added
                 to cutoff lines.}
  \item{text.abline.x}{text to be added to cutoff lines in x direction; if \code{NULL}
   (default) we use ``[pp] \%-cutoff = [ff]'' where [pp] is the percentage up to 2 digits
          and [ff] is the cutoff value up to 2 digits.}
  \item{text.abline.y}{text to be added to cutoff lines in y direction; if \code{NULL}
   (default) we use ``[pp] \%-cutoff = [ff]'' where [pp] is the percentage up to 2 digits
          and [ff] is the cutoff value up to 2 digits.}
  \item{cex.abline}{vector of numerics (cast to length 2): cex-value for added cutoff text.}
  \item{col.abline}{vector of length 2: color for added cutoff text.}
  \item{font.abline}{vector of length 2: font for added cutoff text.}
  \item{adj.abline}{cast to 2 x 2 matrix (by recycling rules): adjustment values for added cutoff text.}
  \item{text.abline.x.y}{y-coordinate of text to be added to cutoff lines in x direction;
                         if \code{NULL} (default) set to mid of \code{mean(par("usr")[c(3,4)])}.}
  \item{text.abline.y.x}{x-coordinate of text to be added to cutoff lines in y direction;
                         if \code{NULL} (default) set to mid of \code{mean(par("usr")[c(1,2)])}.}
  \item{text.abline.x.x}{x-coordinate of text to be added to cutoff lines in x direction;
                         if \code{NULL} (default) set to 1.05 times the cutoff value.}
  \item{text.abline.y.y}{y-coordinate of text to be added to cutoff lines in y direction;
                         if \code{NULL} (default) set to 1.05 times the cutoff value.}
  \item{text.abline.x.fmt.cx}{format string (see \code{\link[base:sprintf]{gettextf}})
        to format the cutoff value in label in x direction.}
  \item{text.abline.x.fmt.qx}{format string to format cutoff probability in label in x direction.}
  \item{text.abline.y.fmt.cy}{format string to format the cutoff value in label in y direction.}
  \item{text.abline.y.fmt.qy}{format string to format cutoff probability in label in y direction.}
  \item{robCov.x}{shall x-distances be based on MCD, i.e.,
                   robust covariances (TRUE) or on classical covariance be used?}
  \item{robCov.y}{shall y-distances be based on MCD, i.e.,
                   robust covariances (TRUE) or on classical covariance be used?}
  \item{tf.x}{transformation for x axis: a function returning the
                    transformed x-coordinates when applied to the data;
              if \code{tf.x} is \code{NULL} (default), internally this
              is set to the evaluation function of the \code{IC.x}.}
  \item{tf.y}{transformation for y axis: a function returning the
          transformed y-coordinates when applied to the data;
          if \code{tf.x} is \code{NULL} (default), internally this
              is set to the evaluation function of \code{IC.y}.}
  \item{jitter.fac}{factor for jittering, see \code{\link{jitter}};}
  \item{jitter.tol}{threshold for jittering: if distance between points is smaller
               than \code{jitter.tol}, points are considered replicates.}
 \item{doplot}{logical; shall a plot be produced? if \code{FALSE} only the return values are produced.}
  \item{main}{the main title.}
}
\details{
calls a corresponding \code{\link{ddPlot}} method to produce the plot.
}
\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }

\value{
If argument \code{doplot} is \code{FALSE}:
A list (returned as \code{invisible()}) with items
\item{id.x}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{x}-cutoff}
\item{id.y}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{y}-cutoff}
\item{id.xy}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{x}-cutoff and the \code{y}-cutoff}
\item{qtx}{the quantiles of the distances of the (possibly transformed) data in \code{x} direction}
\item{qty}{the quantiles of the distances of the (possibly transformed) data in \code{y} direction}
\item{cutoff.x.v}{the cutoff value in \code{x} direction}
\item{cutoff.y.v}{the cutoff value in \code{y} direction}

If argument \code{doplot} is \code{TRUE}:
An S3 object of class \code{c("plotInfo","DiagnInfo")}, i.e., a list
  containing the information needed to produce the
  respective plot, which at a later stage could be used by different
  graphic engines (like, e.g. \code{ggplot}) to produce the plot
  in a different framework. A more detailed description will follow in
  a subsequent version.a list (returned as \code{invisible()}) with items;
  one item is \code{retV} which is the return value in case \code{doplot} is \code{FALSE}.

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

