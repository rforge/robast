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
\S4method{infoPlot}{IC}(object, data = NULL,
             ..., withSweave = getdistrOption("withSweave"),
             col = par("col"), lwd = par("lwd"), lty,
             colI = grey(0.5), lwdI = 0.7*par("lwd"), ltyI = "dotted",
             main = FALSE, inner = TRUE, sub = FALSE,
             col.inner = par("col.main"), cex.inner = 0.8,
             bmar = par("mar")[1], tmar = par("mar")[3],
             with.automatic.grid = TRUE,
             with.legend = TRUE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             x.vec = NULL, scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
             pch.pts = 19,
             cex.npts = 1, cex.npts.fun = NULL, col.npts = grey(.5),
             pch.npts = 20,
             jitter.fac = 1, with.lab = FALSE, cex.lbs = 1, adj.lbs = c(0, 0),
             col.lbs = col.pts, lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
             which.lbs = NULL, which.Order  = NULL, which.nonlbs = NULL,
             attr.pre = FALSE, return.Order = FALSE,
             ylab.abs = "absolute information",
             ylab.rel= "relative information",
             withSubst = TRUE)
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
  \item{tmar}{top margin -- useful for non-standard main title sizes;
              may be a vector with individual values for
              each of the panels to be plotted. }
  \item{bmar}{bottom margin -- useful for non-standard sub title sizes;
              may be a vector with individual values for
              each of the panels to be plotted. }
  \item{col}{color of IC in argument \code{object}.}
  \item{lwd}{linewidth of IC in argument \code{object}.}
  \item{lty}{line-type of IC in argument \code{object}.}
  \item{colI}{color of the classically optimal IC.}
  \item{lwdI}{linewidth of the classically optimal IC.}
  \item{ltyI}{line-type of the classically optimal IC.}
  \item{cex.inner}{magnification to be used for inner titles relative
          to the current setting of \code{cex}; as in 
          \code{\link[graphics]{par}}.}
  \item{col.inner}{character or integer code; color for the inner title}              
  \item{with.automatic.grid}{logical; should a grid be plotted alongside
      with the ticks of the axes, automatically? If \code{TRUE} a respective
      call to \code{grid} in argument \code{panel.first} is ignored. }
  \item{with.legend}{logical; shall a legend be plotted?}
  \item{legend}{either \code{NULL} or a list of length (number of plotted panels)
                of items which can be used as argument \code{legend} in
                command \code{legend}.}
  \item{legend.location}{a valid argument \code{x} for \code{\link{legend}} ---
                         the place where to put the legend on the last issued
                         plot --- or a list of length (number of plotted panels) 
                         of such arguments, one for each plotted panel.}
  \item{legend.bg}{background color for the legend}
  \item{legend.cex}{magnification factor for the legend}
  \item{x.vec}{a numeric vector of grid points to evaluate the influence curve;
               by default, \code{x.vec} is \code{NULL}; then the grid is
               produced automatically according to the distribution of the IC.
               \code{x.vec} can be useful for usage with a rescaling of the
               x-axis to avoid that the evaluation points be selected too
               unevenly (i.e. on an equally spaced grid in the original scale,
               but then, after rescaling non-equally).
               The grid has to be specified in original scale; i.e.; when used
               with rescaling, it should be chosen non-equally spaced. }
  \item{scaleX}{logical; shall X-axis be rescaled (by default according to the cdf of
          the underlying distribution)?}
  \item{scaleY}{logical; shall Y-axis be rescaled for abs.info-plot
                (by default according to a probit scale)?}
  \item{scaleX.fct}{an isotone, vectorized function mapping the domain of the IC
            to [0,1]; if \code{scaleX} is \code{TRUE} and \code{scaleX.fct} is
            missing, the cdf of the underlying observation distribution.}
  \item{scaleX.inv}{the inverse function to \code{scale.fct}, i.e., an isotone,
            vectorized function mapping [0,1] to the domain of the IC
            such that for any \code{x} in the domain,
            \code{scaleX.inv(scaleX.fct(x))==x}; if \code{scaleX} is \code{TRUE}
            and \code{scaleX.inv} is
            missing, the quantile function of the underlying observation distribution. }
  \item{scaleY.fct}{an isotone, vectorized function mapping the
            range of the norm of the IC to [0,1]; defaulting
            to the cdf of \eqn{{\cal N}(0,1)}{N(0,1)};
            can also be a list of functions with one list element for each
            of the panels to be plot.}
  \item{scaleY.inv}{an isotone, vectorized function mapping [0,1] into the range
            of the norm of the IC; defaulting to the quantile function
            of  \eqn{{\cal N}(0,1)}{N(0,1)};
            can also be a list of functions with one list element for each
            of the panels to be plot.}
  \item{scaleN}{integer; defaults to 9; on rescaled axes, number of x
                and y ticks if drawn automatically;}
  \item{x.ticks}{numeric; defaults to NULL; (then ticks are chosen automatically);
                 if non-NULL, user-given x-ticks (on original scale);}
  \item{y.ticks}{numeric; defaults to NULL; (then ticks are chosen automatically);
                 if non-NULL, user-given y-ticks (on original scale);
                 can be a list with one (numeric or NULL) item per panel}
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
  \item{withSubst}{logical; if \code{TRUE} (default) pattern substitution for
        titles and lables is used; otherwise no substitution is used. }
  \item{col.pts}{color of the points of the \code{data} argument plotted;
     can be a vector or a matrix. More specifically, if argument \code{attr.pre}
     is \code{TRUE}, it is recycled to fill a matrix of dimension \code{n} by 2
     (\code{n} the number of observations prior to any selection) where filling
     is done in order column first. The two columns are used for possibly
     different colors for the actual IC from the argument and the classical
     IC which is also shown. The selection done via \code{which.lbs} and
     \code{which.Order} is then done afterwards and on this matrix;
     argument \code{col.npts} is ignored in this case. If \code{attr.pre} is \code{FALSE},
     \code{col.pts} is recycled to fill a matrix of dimension \code{n.s} by 2
     where \code{n.s} is the number of observations selected for labelling
     and refers to the index ordering after the selection. Then argument
     \code{col.npts} deteremines the colors of the shown but non-labelled
     observations as given in argument \code{which.nonlbs}.}
  \item{pch.pts}{symbol of the points of the \code{data} argument plotted
                 (may be a vector of length 2 or a matrix, see \code{col.pts},
                 with argument \code{pch.npts} as counterpart).}
  \item{cex.pts}{size of the points of the \code{data} argument plotted
                 (may be a vector of length 2 or a matrix, see \code{col.pts},
                 with argument \code{cex.npts} as counterpart).}
  \item{cex.pts.fun}{rescaling function for the size of the points to be plotted;
        either \code{NULL} (default), then \code{log(1+abs(x))} is used for each of
        the rescalings, or a function which is then used for each of the
        rescalings, or a list of functions; if it is a function or a list of
        functions, if necessary it is recylced to length \code{2 * dim}
        where \code{2} is for the classical IC and the IC in argument \code{object}
        and \code{dim} is the number of dimensions of the pICs to be plotted;
        in the index of this list, \code{2} is incremented first;
        then \code{dim}.}
  \item{col.npts}{color of the non-labelled points of the \code{data} argument
     plotted; (may be a vector of length 2, or it can be a matrix
     \code{nnlb <- sum(which.nonlbs)} by \code{2},
     \code{nnlb} the number of non-labelled shown observations.}
  \item{pch.npts}{symbol of the non-labelled points of the \code{data} argument
         plotted (may be a vector of length \code{2} or a matrix,
         see \code{col.npts}).}
  \item{cex.npts}{size of the non-labelled points of the \code{data} argument
                 plotted (may be a vector of length \code{2} or a matrix, see
                 \code{col.npts}).}
  \item{cex.npts.fun}{rescaling function for the size of the non-labelled points
  to be plotted; either \code{NULL} (default), then \code{log(1+abs(x))}
  is used for each of the rescalings, or a function which is then used
  for each of the rescalings, or a list of functions; if it is a
  function or a list of functions, if necessary it is recylced
  to length \code{2 * dim} where \code{dim} is the number of dimensions of
  the pICs to be plotted;  in the index of this list,
  \code{2} is incremented first; then \code{dim}.}
  \item{attr.pre}{logical; do graphical attributes for plotted data refer
                  to indices prior (\code{TRUE}) or posterior to selection
                  via arguments \code{which.lbs}, \code{which.Order}, \code{which.nonlbs}
                  (\code{FALSE})? }
  \item{with.lab}{logical; shall labels be plotted to the observations?
                  (may be a vector of length 2, see \code{col.pts} -- 
                   but not a matrix)}
  \item{cex.lbs}{size of the labels; can be vectorized to an array
                 of dim nlbs x 2 x npnl where npnl is the number of plotted
                 panels and nlbs the number of plotted labels; if it is
                 a vector, it is recylced in order labels then ICs
                 [arg \code{IC}/classic] then panels.}
  \item{col.lbs}{color of the labels; can be vectorized to a matrix
                 of dim nlbs x 2 as \code{col.pts}.}
  \item{adj.lbs}{adjustment of the labels; can be vectorized to an array
             of dim 2 x 2 x npnl matrix, npnl the number of plotted panels;
             if it is a vector, it is recycled in order (x,y)-coords then ICs
                 [arg \code{IC}/classic] then panels.}
  \item{lab.pts}{character or NULL; labels to be plotted to the observations;
                 can be a vector of length \code{n}, \code{n} the number of 
                 all observations prior to any selection with \code{which.lbs},
                 \code{which.Order}; if \code{lab.pts} is \code{NULL}, 
                 observation indices are used.}
  \item{lab.font}{font to be used for labels; (may be a vector of length 2, 
                  see \code{with.lab}).}
  \item{alpha.trsp}{alpha transparency to be added ex post to colors
        \code{col.pch} and \code{col.nonlbl}; if one-dim and NA all colors are
        left unchanged. Otherwise, with usual recycling rules \code{alpha.trsp}
        gets shorted/prolongated to length the number of panel data-symbols to 
        be plotted. Coordinates of this vector \code{alpha.trsp} with NA are left unchanged,
        while for the remaining ones, the alpha channel in rgb space is set
        to the respective coordinate value of \code{alpha.trsp}. The non-NA
        entries must be integers in [0,255] (0 invisible, 255 opaque).}
  \item{jitter.fac}{jittering factor used in case of a \code{DiscreteDistribution}
                    for plotting points of the \code{data} argument in a jittered fashion
                  (may be a vector of length 2, see \code{with.lab}).}
  \item{which.lbs}{either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then no observation is excluded}
  \item{which.Order}{we order the observations (descending) according to the norm given by
           \code{normtype(object)}; then \code{which.Order}
           either is an integer vector with the indices of the \emph{ordered}
           observations (remaining after a possible reduction by argument \code{which.lbs})
           to be plotted into graph or \code{NULL} --- then no (further) observation
           is excluded.}
  \item{which.nonlbs}{indices of the observations which should be plotted but
        not labelled; either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then all non-labelled
          observations are plotted.}
  \item{return.Order}{logical; if \code{TRUE}, a list of length two with order vectors
    is returned --- one for ordering w.r.t. the given IC, one for ordering w.r.t.
    the classically optimal IC; more specifically, the order of the (remaining) observations
    given by their original index is returned (remaining means: after a possible
    reduction by argument \code{which.lbs}, and ordering is according to the norm given by
           \code{normtype(object)});
   otherwise we return \code{invisible()} as usual.}
  \item{ylab.abs}{character; label to be used for y-axis in absolute information panel} 
  \item{ylab.rel}{character; label to be used for y-axis in relative information panel} 
  \item{\dots}{further parameters for \code{plot}}
}
\value{An S3 object of class \code{c("plotInfo","DiagnInfo")}, i.e., a list
  containing the information needed to produce the
  respective plot, which at a later stage could be used by different
  graphic engines (like, e.g. \code{ggplot}) to produce the plot
  in a different framework. A more detailed description will follow in
  a subsequent version.
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
has correct length. If argument \code{withSubst} is \code{TRUE}, in all title 
and axis lable arguments, the following patterns are substituted:
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

Similarly, if argument \code{\dots} contains arguments \code{xaxt} or
\code{yaxt}, these may be vectorized, with one value for each of the panels
to be plotted. This is useful for stacking panels over each other, using
a common x-axis (see example below).

The \code{\dots} argument may also contain an argument \code{withbox} which
if \code{TRUE} warrants that even if \code{xaxt} and \code{yaxt} both are
\code{FALSE}, a box is drawn around the respective panel.

In addition, argument \code{\dots} may contain arguments \code{panel.first},
\code{panel.last}, i.e., hook expressions to be evaluated at the very beginning
and at the very end of each panel (within the then valid coordinates).
To be able to use these hooks for each panel individually, they may also be
lists of expressions (of the same length as the number of panels and
run through in the same order as the panels).
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

## don't run to reduce check time on CRAN
\donttest{
## selection of subpanels for plotting
par(mfrow=c(1,2))
infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","sd"))
infoPlot(IC1, mfColRow = FALSE, to.draw.arg=c("Abs","sd"), log="y")

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

ICr <- makeIC(list(function(x)sign(x),function(x)sign(abs(x)-qnorm(.75))),N)
data <- r(N)(600)
data.c <- c(data, 1000*data[1:30])
par(mfrow=c(3,1))
infoPlot(ICr, data=data.c, tmar=c(4.1,0,0), bmar=c(0,0,4.1),
         xaxt=c("n","n","s"), mfColRow = FALSE, panel.first= grid(),
         cex.pts=c(.9,.9), alpha.trsp=20, lwd=2, lwdI=1.5, col=3,
         col.pts=c(3,2), colI=2, pch.pts=c(20,20), inner=FALSE,
         scaleX = TRUE, scaleX.fct=pnorm, scaleX.inv=qnorm,
         scaleY=TRUE, scaleY.fct=function(x) pchisq(x,df=1),
         scaleY.inv=function(x)qchisq(x,df=1),legend.cex = 1.0)

}

}
\concept{absolute information}
\concept{relative information}
\keyword{robust}
