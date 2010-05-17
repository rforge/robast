\name{infoPlot}
\alias{infoPlot}
\alias{infoPlot-methods}
\alias{infoPlot,IC-method}

\title{Plot absolute and relative information}
\description{
  Plot absolute and relative information of influence curves.
}
\usage{
infoPlot(object,  ...)
\S4method{infoPlot}{IC}(object, data = NULL, ...,
             withSweave = getdistrOption("withSweave"),
             col = par("col"), lwd = par("lwd"), lty, 
             colI = grey(0.5), lwdI = 0.7*par("lwd"), ltyI = "dotted",
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3],
             with.legend = TRUE, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, col.pts = par("col"),
             pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
             lab.pts = NULL, lab.font = NULL,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE,
             ylab.abs = "absolute information", 
             ylab.rel= "relative information")
}
\arguments{
  \item{object}{object of class \code{"InfluenceCurve"} }
  \item{data}{optional data argument --- for plotting observations into the plot;}
  \item{withSweave}{logical: if \code{TRUE} (for working with \command{Sweave}) 
        no extra device is opened}
  \item{main}{logical: is a main title to be used? or \cr
              just as argument \code{main} in \code{\link{plot.default}}.}
  \item{inner}{logical: do panels have their own titles? or \cr
               character vector of / cast to length 'number of compared 
               dimensions';
               if argument \code{to.draw.arg} is used, this refers to 
               a vector of length \code{1} (absolute information) 
               \code{+ length(to.draw.arg)}, the actually plotted
               relative informations. For further information, see  also
               \code{main} in \code{\link{plot.default}}.} 
  \item{sub}{logical: is a sub-title to be used? or \cr
              just as argument \code{sub} in \code{\link{plot.default}}.}
  \item{tmar}{top margin -- useful for non-standard main title sizes}
  \item{bmar}{bottom margin -- useful for non-standard sub title sizes}
  \item{col}{color of IC in argument \code{object}.}
  \item{lwd}{linewidth of IC in argument \code{object}.}
  \item{lty}{line-type of IC in argument \code{object}.}
  \item{colI}{color of the classically optimal IC}
  \item{lwdI}{linewidth of the classically optimal IC}
  \item{ltyI}{line-type of the classically optimal IC}
  \item{cex.inner}{magnification to be used for inner titles relative
          to the current setting of \code{cex}; as in 
          \code{\link[graphics]{par}}}
  \item{col.inner}{character or integer code; color for the inner title}              
  \item{with.legend}{logical; shall a legend be plotted?}
  \item{legend.location}{a valid argument \code{x} for \code{\link{legend}} ---
                         the place where to put the legend on the last issued
                         plot --- or a list of length (number of plotted panels) 
                         of such arguments, one for each plotted panel.}
  \item{legend.bg}{background color for the legend}
  \item{legend.cex}{magnification factor for the legend}
  \item{mfColRow}{shall default partition in panels be used --- defaults to \code{TRUE}}
  \item{to.draw.arg}{Either \code{NULL} (default; 
                          everything is plotted) or a vector making a selection
                          among the relative information plots; the absolute
                          information being plotted in any case. This
                          vector is either a vector of integers 
                         (the indices of the subplots to be drawn) or characters 
                         --- the names of the subplots to be drawn: these
                         names are to be chosen either among the row names of 
                         the trafo matrix 
                         \code{rownames(trafo(eval(object@CallL2Fam)@param))} 
                         or if the last expression is \code{NULL} a 
                         vector \code{"dim<dimnr>"}, \code{dimnr} running through 
                         the number of rows of the trafo matrix.
                         }
  \item{cex.pts}{size of the points of the \code{data} argument plotted}
  \item{col.pts}{color of the points of the \code{data} argument plotted}
  \item{pch.pts}{symbol of the points of the \code{data} argument plotted}
  \item{with.lab}{logical; shall labels be plotted to the observations?}
  \item{lab.pts}{character or NULL; labels to be plotted to the observations; if \code{NULL}
                 observation indices;}
  \item{lab.font}{font to be used for labels}
  \item{jitter.fac}{jittering factor used in case of a \code{DiscreteDistribution}
                    for plotting points of the \code{data} argument in a jittered fashion.}
  \item{which.lbs}{either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then no observation is excluded}
  \item{which.Order}{we order the observations (descending) according to the norm given by
           \code{normtype(object)}; then \code{which.Order}
           either is an integer vector with the indices of the \emph{ordered}
           observations (remaining after a possible reduction by argument \code{which.lbs})
           to be plotted into graph or \code{NULL} --- then no (further) observation
           is excluded.}
  \item{return.Order}{logical; if \code{TRUE}, a list of length two with order vectors
    is returned --- one for ordering w.r.t. the given IC, one for ordering w.r.t.
    the classically optimal IC; more specifically, the order of the (remaining) observations
    given by their original index is returned (remaining means: after a possible
    reduction by argument \code{which.lbs}, and ordering is according to the norm given by
           \code{normtype(object)});
   othervise we return \code{invisible()} as usual.}
  \item{ylab.abs}{character; label to be used for y-axis in absolute information panel} 
  \item{ylab.rel}{character; label to be used for y-axis in relative information panel} 
  \item{\dots}{further parameters for \code{plot}}
}
\details{
  Absolute information is defined as the square of the length
  of an IC. The relative information is defined as the 
  absolute information of one component with respect to the 
  absolute information of the whole IC; confer Section 8.1 
  of Kohl (2005). \cr

Any parameters of \code{plot.default} may be passed on to this particular
\code{plot} method. 

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
length 2*(number of plotted dimensions + e), where e is 1 or 0 depending
on whether absolute information is plotted or not; 
in the case of longer length, 
if e is 1, the first two elements are the values for \code{ylim} 
in panel "Abs", while the last 2*(number of plotted dimensions)
are the values for \code{ylim} for the plotted dimensions of the IC, 
one pair for each dimension.
}
%\value{}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}}
\examples{
N <- NormLocationScaleFamily(mean=0, sd=1) 
IC1 <- optIC(model = N, risk = asCov())
infoPlot(IC1)
## selection of subpanels for plotting
par(mfrow=c(1,2))
infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","sd"))
infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","mean"), 
              panel.first= grid(), ylim = c(0,4), xlim = c(-6,6))
infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","mean"), 
              panel.first= grid(), ylim = c(0,4,-3,3), xlim = c(-6,6))
par(mfrow=c(1,3))
infoPlot(IC1, mfColRow = FALSE, panel.first= grid(),
         ylim = c(0,4,0,.3,0,.8), xlim=c(-6,6))
par(mfrow=c(1,1))
data <- r(N)(20)
par(mfrow=c(1,3))
infoPlot(IC1, data=data, mfColRow = FALSE, panel.first= grid(),
         with.lab = TRUE, cex.pts=2,
         which.lbs = c(1:4,15:20), which.Order = 1:6,
         return.Order = TRUE)
infoPlot(IC1, data=data[1:10], mfColRow = FALSE, panel.first= grid(),
         with.lab = TRUE, cex.pts=0.7)
par(mfrow=c(1,1))
}
\concept{absolute information}
\concept{relative information}
\keyword{robust}
