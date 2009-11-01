\name{internals_for_RobAStBase_ddPlot}
\alias{internals_for_RobAStBase_ddPlot}
\alias{.ddPlot.MatNtNtCoCo}

\title{Internal / Helper functions of package RobAStBase for ddPlot}

\description{
This function is an internally used helper function for \code{\link{ddPlot}} in package \pkg{RobAStBase}.}

\usage{
.ddPlot.MatNtNtCoCo(data, ...,  dist.x = NormType(), dist.y  = NormType(),
       cutoff.x = cutoff(norm = dist.x, cutoff.quantile  = cutoff.quantile.x),
       cutoff.y = cutoff(norm = dist.y, cutoff.quantile  = cutoff.quantile.y),
       cutoff.quantile.x = 0.95, cutoff.quantile.y = cutoff.quantile.x,
       transform.x, transform.y = transform.x,
       id.n, lab.pts, adj =0, cex.idn = 1,
       col.idn = par("col"), lty.cutoff,
       lwd.cutoff, col.cutoff = "red")}

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
