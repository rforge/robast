\name{internals_for_RobAStBase_ddPlot}
\alias{internals_for_RobAStBase_ddPlot}
\alias{.ddPlot.MatNtNtCoCo}

\title{Internal / Helper functions of package RobAStBase for ddPlot}

\description{
This function is an internally used helper function for \code{\link{ddPlot}} in package \pkg{RobAStBase}.}

\usage{
.ddPlot.MatNtNtCoCo(data, ..., dist.x = NormType(), dist.y  = NormType(),
cutoff.x = cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x),
cutoff.y = cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y),
cutoff.quantile.x = 0.95,  cutoff.quantile.y = cutoff.quantile.x,
transform.x,transform.y = transform.x, id.n, cex.pts = 1, lab.pts,
jitt.pts = 0, alpha.trsp = NA, adj =0,cex.idn = 1,col.idn = par("col"),
lty.cutoff,lwd.cutoff,col.cutoff = "red",text.abline = TRUE,
text.abline.x = NULL, text.abline.y = NULL,cex.abline = par("cex"), 
col.abline = col.cutoff,font.abline = par("font"), adj.abline = c(0,0),
text.abline.x.x = NULL, text.abline.x.y = NULL, text.abline.y.x = NULL, 
text.abline.y.y = NULL, text.abline.x.fmt.cx = "\%7.2f",
text.abline.x.fmt.qx = "\%4.2f\%\%", text.abline.y.fmt.cy = "\%7.2f",
text.abline.y.fmt.qy = "\%4.2f\%\%", jitt.fac = 10, doplot = TRUE)
}
\arguments{
  \item{data}{data in \code{matrix} form (columns are observations; rows are variable
  dimensions) at which to produce the \code{ddPlot}.}
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
  \item{cex.pts}{the corresponding \code{cex} argument for plotted points.}
  \item{lab.pts}{a vector of labels for the (unsubsetted) \code{data}.}
  \item{jit.pts}{the corresponding \code{jitter} argument for plotted points;
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
  \item{text.abline.x.fmt.cx}{format string (see \code{\link[base:sprintf]{gettextf}}) to format the cutoff value in label in x direction.}
  \item{text.abline.x.fmt.qx}{format string to format cutoff probability in label in x direction.}
  \item{text.abline.y.fmt.cy}{format string to format the cutoff value in label in y direction.}
  \item{text.abline.y.fmt.qy}{format string to format cutoff probability in label in y direction.}
  \item{jit.fac}{factor for jittering, see \code{jitter};}
  \item{jit.tol}{threshold for jittering: if distance between points is smaller
               than \code{jit.tol}, points are considered replicates.}
  \item{doplot}{logical; shall a plot be produced? if \code{FALSE} only the return values are produced.}
}

\details{
performs the plotting for \code{\link{ddPlot}} and \code{\link{outlyingPlotIC}};
all arguments except for \code{data} are optional. In case they are missing
default values are used as usual; for those arguments without default arguments,
we do
\describe{
\item{transform.x}{defaults to identity, internally}.
\item{id.n}{defaults to \code{1:ncol(data)}, internally}.
\item{lab.pts}{defaults to \code{(1:ncol(data))[id.n]}, internally}.
\item{lwd.cutoff}{defaults to argument \code{lwd}, if given, else to \code{par{lwd}}, internally}.
\item{lty.cutoff}{defaults to argument \code{lty}, if given, else to \code{par{lty}}, internally}.
}
}


\value{
a list (returned as \code{invisible()}) with items
\item{id.x}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{x}-cutoff}
\item{id.y}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{y}-cutoff}
\item{id.xy}{the indices of (possibly transformed) data (within subset \code{id.n}) beyond the \code{x}-cutoff and the \code{y}-cutoff}
\item{qtx}{the quantiles of the distances of the (possibly transformed) data in \code{x} direction}
\item{qty}{the quantiles of the distances of the (possibly transformed) data in \code{y} direction}
\item{cutoff.x.v}{the cutoff value in \code{x} direction}
\item{cutoff.y.v}{the cutoff value in \code{y} direction}
}


\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  }
\seealso{\code{\link[graphics]{plot.default}}, \code{\link[graphics]{par}},
          \code{\link{ddPlot}}, \code{\link{outlyingPlotIC}}}

\examples{
MX <- matrix(rnorm(1500),nrow=6)
QM <- matrix(rnorm(36),nrow=6); QM <- QM \%*\% t(QM)
RobAStBase:::.ddPlot.MatNtNtCoCo(data=MX, 
        dist.y=QFNorm(QuadF=PosSemDefSymmMatrix(QM)),
        xlab="Norm.x",ylab="Norm.y", cex.idn = 1.3, offset=0,
        lwd=2, lwd.cutoff=4, lty=2, col.cutoff =2, col.idn="green",
        col = "blue", adj=0.4, pos=4,id.n = sample(1:200,size=100),
        lab.pts=letters,log="x", main="GA", sub="NO",cex.sub=0.2)
}
\keyword{internal}
\concept{utilities}
\keyword{hplot}
