\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,IC,missing-method}
\alias{plot,IC,numeric-method}
\title{ Methods for Function plot in Package `RobAStBase' }
\description{plot-methods}
\usage{plot(x, y, ...)
\S4method{plot}{IC,missing}(x, ..., withSweave = getdistrOption("withSweave"),
             main = FALSE, inner = TRUE, sub = FALSE,
             col.inner = par("col.main"), cex.inner = 0.8,
             bmar = par("mar")[1], tmar = par("mar")[3],
             mfColRow = TRUE, to.draw.arg = NULL)
\S4method{plot}{IC,numeric}(x, y, ..., cex.pts = 1,
                            col.pts = par("col"), pch.pts = 1, jitter.fac = 1,
                            with.lab = FALSE, lab.pts = NULL, lab.font = NULL,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE)
}
\arguments{
  \item{x}{object of class \code{"IC"}: IC to be plotted }
  \item{y}{ missing or numeric (a dataset, e.g.)}
  \item{withSweave}{logical: if \code{TRUE} (for working with \command{Sweave})
        no extra device is opened}
  \item{main}{logical: is a main title to be used? or \cr
              just as argument \code{main} in \code{\link{plot.default}}.}
  \item{inner}{logical: do panels have their own titles? or \cr
               character vector of / cast to length 'number of plotted 
               dimensions';
               if argument \code{to.draw.arg} is used, this refers to 
               a vector of length \code{length(to.draw.arg)}, the 
               actually plotted dimensions. For further information, see also 
               description of argument \code{main} in \code{\link{plot.default}}.} 
  \item{sub}{logical: is a sub-title to be used? or \cr
              just as argument \code{sub} in \code{\link{plot.default}}.}
  \item{tmar}{top margin -- useful for non-standard main title sizes}
  \item{bmar}{bottom margin -- useful for non-standard sub title sizes}
  \item{cex.inner}{magnification to be used for inner titles relative
          to the current setting of \code{cex}; as in
          \code{\link[graphics]{par}}}
  \item{col.inner}{character or integer code; color for the inner title}
  \item{mfColRow}{shall default partition in panels be used --- defaults to \code{TRUE}}
  \item{to.draw.arg}{Either \code{NULL} (default; 
                          everything is plotted) or a vector of either integers 
                         (the indices of the subplots to be drawn) or characters 
                         --- the names of the subplots to be drawn: these
                         names are to be chosen either among the row names of 
                         the trafo matrix 
                         \code{rownames(trafo(eval(x@CallL2Fam)@param))} 
                         or if the last expression is \code{NULL} a 
                         vector \code{"dim<dimnr>"}, \code{dimnr} running through 
                         the number of rows of the trafo matrix.
                         }
  \item{cex.pts}{size of the points of the second argument plotted}
  \item{col.pts}{color of the points of the second argument plotted}
  \item{pch.pts}{symbol of the points of the second argument plotted}
  \item{with.lab}{logical; shall labels be plotted to the observations?}
  \item{lab.pts}{character or NULL; labels to be plotted to the observations; if \code{NULL}
                 observation indices;}
  \item{lab.font}{font to be used for labels}
  \item{jitter.fac}{jittering factor used in case of a \code{DiscreteDistribution}
                    for plotting points of the second argument in a jittered fashion.}
  \item{which.lbs}{either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then no observation is excluded}
  \item{which.Order}{we order the observations (descending) according to the norm given by
           \code{normtype(object)}; then \code{which.Order}
           either is an integer vector with the indices of the \emph{ordered}
           observations (remaining after a possible reduction by argument \code{which.lbs})
           to be plotted into graph or \code{NULL} --- then no (further) observation
           is excluded.}
  \item{return.Order}{logical; if \code{TRUE}, an order vector
    is returned; more specifically, the order of the (remaining) observations
    given by their original index is returned (remaining means: after a possible
    reduction by argument \code{which.lbs}, and ordering is according to the norm given by
           \code{normtype(object)});
   othervise we return \code{invisible()} as usual.}
  \item{\dots}{further parameters for \code{plot}}
}
\details{
Any parameters of \code{plot.default} may be passed on to this particular
\code{plot} method.

We start describing the \code{IC,missing}-method:
For main-, inner, and subtitles given as arguments \code{main},
\code{inner}, and \code{sub}, top and bottom margins are enlarged to 5 resp.
6 by default but may also be specified by \code{tmar} / \code{bmar} arguments.
If \code{main} / \code{inner} / \code{sub} are
logical then if the respective argument is \code{FALSE} nothing is done/plotted,
but if it is \code{TRUE}, we use a default main title taking up the calling
arguments in case of \code{main}, default inner titles taking up the
class and (named) parameter slots of arguments in case of \code{inner},
and a "generated on <data>"-tag in case of \code{sub}.
Of course, if \code{main} / \code{inner} / \code{sub} are \code{character}, this
is used for the title; in case of \code{inner} it is then checked whether it
has correct length. In all title arguments, the following patterns are substituted:
\describe{
\item{\code{"\%C"}}{class of argument \code{object}}
\item{\code{"\%A"}}{deparsed argument  \code{object}}
\item{\code{"\%D"}}{time/date-string when the plot was generated}
}

If argument \code{\dots} contains argument \code{ylim}, this may either be
as in \code{plot.default} (i.e. a vector of length 2) or a vector of 
length 2*(number of plotted dimensions + 2), 
where the first two elements are the values for \code{ylim} 
in panel "d",  the first two are for \code{ylim} resp. \code{xlim} 
for panels "p" and "q", and the last 2*(number of plotted dimensions)
are the values for \code{ylim} for the plotted dimensions of the L2derivative, 
one pair for each dimension.

The \code{IC,numeric}-method calls the \code{IC,missing}-method but in
addition plots the values of a dataset into the IC.

}
\examples{
IC1 <- new("IC")
plot(IC1)
plot(IC1, main = TRUE, panel.first= grid(),
     col = "blue", cex.main = 2, cex.inner = 1)

### selection of subpanels for plotting
N <- NormLocationScaleFamily(mean=0, sd=1) 
IC2 <- optIC(model = N, risk = asCov())
par(mfrow=c(1,1))
plot(IC2, main = TRUE, panel.first= grid(),
     col = "blue", cex.main = 2, cex.inner = 0.6,
     mfColRow = FALSE, to.draw.arg=c("sd"))

## xlim and ylim arguments
plot(IC2, main = TRUE, panel.first= grid(), 
     ylim=c(-3,3), xlim=c(-2,3))
plot(IC2, main = TRUE, panel.first= grid(), 
     ylim=c(-3,3,-1,3), xlim=c(-2,3))

data <- r(N)(30)
plot(IC2, data, panel.first= grid(),
     ylim = c(-3,3,-1,3), xlim=c(-2,3),
     cex.pts = 3, pch.pts = 1:2, col.pts="green",
     with.lab = TRUE, which.lbs = c(1:4,15:20),
     which.Order = 1:6, return.Order = TRUE)
}
\keyword{methods}
\keyword{distribution}

