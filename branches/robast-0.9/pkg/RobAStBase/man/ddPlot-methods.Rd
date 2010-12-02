\name{ddPlot-methods}
\docType{methods}
\alias{ddPlot}
\alias{ddPlot-methods}
\alias{ddPlot,matrix-method}
\alias{ddPlot,numeric-method}
\alias{ddPlot,data.frame-method}
\title{ Methods for Function ddPlot in Package `RobAStBase' }
\description{ddPlot-methods}
\usage{
ddPlot(data, dist.x, dist.y, cutoff.x, cutoff.y, ...)
\S4method{ddPlot}{matrix}(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL, 
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "\%7.2f", text.abline.x.fmt.qx = "\%4.2f\%\%",
       text.abline.y.fmt.cy = "\%7.2f", text.abline.y.fmt.qy = "\%4.2f\%\%")
\S4method{ddPlot}{numeric}(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL, 
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "\%7.2f", text.abline.x.fmt.qx = "\%4.2f\%\%",
       text.abline.y.fmt.cy = "\%7.2f", text.abline.y.fmt.qy = "\%4.2f\%\%")
\S4method{ddPlot}{data.frame}(data, dist.x = NormType(), dist.y  = NormType(),
       cutoff.x, cutoff.y, ...,
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj, cex.idn,
       col.idn, lty.cutoff, lwd.cutoff, col.cutoff, text.abline = TRUE,
       text.abline.x = NULL, text.abline.y = NULL,
       cex.abline = par("cex"), col.abline = col.cutoff,
       font.abline = par("font"), adj.abline = c(0,0),
       text.abline.x.x = NULL, text.abline.x.y = NULL, 
       text.abline.y.x = NULL, text.abline.y.y = NULL,
       text.abline.x.fmt.cx = "\%7.2f", text.abline.x.fmt.qx = "\%4.2f\%\%",
       text.abline.y.fmt.cy = "\%7.2f", text.abline.y.fmt.qy = "\%4.2f\%\%")
}
\arguments{
  \item{data}{data coercable to \code{matrix}; the data  at which to produce the \code{ddPlot}.}
  \item{\dots}{further arguments to be passed to \code{plot.default}, \code{text}, and \code{abline}}
  \item{dist.x}{object of class \code{NormType}; the distance for the \code{x} axis.}
  \item{dist.y}{object of class \code{NormType}; the distance for the \code{y} axis.}
  \item{cutoff.x}{object of class \code{cutoff}; the cutoff information for the \code{x} axis
  (the vertical line discriminating 'good' and 'bad' points).}
  \item{cutoff.y}{object of class \code{cutoff}; the cutoff information for the \code{y} axis
  (the horizontal line discriminating 'good' and 'bad' points).}
  \item{cutoff.quantile.x}{numeric; the cutoff quantile for the \code{x} axis.}
  \item{cutoff.quantile.y}{numeric; the cutoff quantile for the \code{y} axis.}
  \item{transform.x}{function; a transformation to be performed before determining  the
  distances of the \code{x} axis.}
  \item{transform.y}{function; a transformation to be performed before determining  the
  distances of the \code{y} axis.}
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
                 \code{\link[graphics]{abline}} for drawing the cutoff lines;
                 either one lty-value (one value or vector) or a list of length 2
                 of lty-values.}
  \item{lwd.cutoff}{(vector cast to length 2): the corresponding \code{lwd} argument for
                 \code{\link[graphics]{abline}} for drawing the cutoff lines.}
  \item{col.cutoff}{(vector cast to length 2): the corresponding \code{col} argument for
                 \code{\link[graphics]{abline}} for drawing the cutoff lines.}
  \item{text.abline}{vector of logicals (cast to length 2): shall text be added
                 to cutoff lines.}
  \item{text.abline.x}{text to be added to cutoff lines in x direction; if \code{NULL}
   (default) we use ``[pp] \%-cutoff = [ff]'' where [pp] is the percentag up to 2 digits
          and [ff] is the cutoff value up to 2 digits.}
  \item{text.abline.y}{text to be added to cutoff lines in y direction; if \code{NULL}
   (default) we use ``[pp] \%-cutoff = [ff]'' where [pp] is the percentag up to 2 digits
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
  \item{text.abline.x.fmt.cx}{format string (see \code{\link[base]{gettextf}}) to format the cutoff value in label in x direction.}
  \item{text.abline.x.fmt.qx}{format string to format cutoff probability in label in x direction.}
  \item{text.abline.y.fmt.cy}{format string to format the cutoff value in label in y direction.}
  \item{text.abline.y.fmt.qy}{format string to format cutoff probability in label in y direction.}
}
\details{
The \code{matrix}-method calls \code{.ddPlot.MatNtNtCoCo},
the \code{numeric}- and \code{data.frame}-methods coerce argument \code{data}
to \code{matrix} --- the \code{numeric}-method by a call to \code{matrix(data, nrow=1)},
in the  \code{data.frame}-methods by a call to \code{t(as.matrix(data))}.

In arguments \code{text.abline.x} and \code{text.abline.y} the following 
patterns are substituted:
\describe{
\item{\code{"\%qx"}}{cutoff-quantile in x-direction}
\item{\code{"\%qy"}}{cutoff-quantile in y-direction}
\item{\code{"\%cx"}}{cutoff-value in x-direction}
\item{\code{"\%cy"}}{cutoff-value in y-direction}
}
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
\author{
  Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}
  }

\examples{
MX <- matrix(rnorm(1500),nrow=6)
QM <- matrix(rnorm(36),nrow=6); QM <- QM \%*\% t(QM)
ddPlot(data=MX, dist.y=QFNorm(QuadF=PosSemDefSymmMatrix(QM)))
}
\keyword{methods}
\keyword{hplot}

