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
             with.automatic.grid = TRUE,
             with.legend = FALSE, legend = NULL, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             withMBR = FALSE, MBRB = NA, MBR.fac = 2, col.MBR = par("col"),
             lty.MBR = "dashed", lwd.MBR = 0.8,
             x.vec = NULL, scaleX = FALSE, scaleX.fct, scaleX.inv,
             scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
             scaleN = 9, x.ticks = NULL, y.ticks = NULL,
             mfColRow = TRUE, to.draw.arg = NULL,
             withSubst = TRUE)
\S4method{plot}{IC,numeric}(x, y, ..., cex.pts = 1, col.pts = par("col"), pch.pts
                 = 1, jitter.fac = 1, with.lab = FALSE, lab.pts = NULL,
                 lab.font = NULL, alpha.trsp = NA, which.lbs = NULL,
                 which.Order = NULL, return.Order = FALSE)
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
  \item{withMBR}{logical; shall horizontal lines with min and max of MBRE be plotted for
                 comparison?}
  \item{MBRB}{matrix (or \code{NA}); coerced by usual recycling rules to a
              matrix with as many rows as plotted panels and with first column
              the lower bounds and the second column the upper bounds for the
              respective coordinates (ideally given by the MBR-IC).}
  \item{MBR.fac}{positive factor; scales the bounds given by argument \code{MBRB}}
  \item{col.MBR}{color for the MBR lines; as usual \code{col}-argument;}
  \item{lty.MBR}{line type for the MBR lines; as usual \code{lty}-argument;}
  \item{lwd.MBR}{line width for the MBR lines; as usual \code{lwd}-argument;}
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
  \item{scaleY}{logical; shall Y-axis be rescaled (by default according to a probit scale)?}
  \item{scaleX.fct}{an isotone, vectorized function mapping the domain of the IC
            to [0,1]; if \code{scaleX} is \code{TRUE} and \code{scaleX.fct} is
            missing, the cdf of the underlying observation distribution;
            can also be a list of functions with one list element for each
            of the panels to be plot.}
  \item{scaleX.inv}{the inverse function to \code{scale.fct}, i.e., an isotone,
            vectorized function mapping [0,1] to the domain of the IC
            such that for any \code{x} in the domain,
            \code{scaleX.inv(scaleX.fct(x))==x}; if \code{scaleX} is \code{TRUE}
            and \code{scaleX.inv} is
            missing, the quantile function of the underlying observation distribution;
            can also be a list of functions with one list element for each
            of the panels to be plot.}
  \item{scaleY.fct}{an isotone, vectorized function mapping for each coordinate the
            range of the respective coordinate of the IC
            to [0,1]; defaulting to the cdf of \eqn{{\cal N}(0,1)}{N(0,1)}.}
  \item{scaleY.inv}{an isotone, vectorized function mapping for each coordinate
            the range [0,1] into the range of the respective coordinate of the IC;
            defaulting to the quantile function of  \eqn{{\cal N}(0,1)}{N(0,1)}.}
  \item{scaleN}{integer; defaults to 9; on rescaled axes, number of x
                and y ticks if drawn automatically;}
  \item{x.ticks}{numeric; defaults to NULL; (then ticks are chosen automatically);
                 if non-NULL, user-given x-ticks (on original scale);}
  \item{y.ticks}{numeric; defaults to NULL; (then ticks are chosen automatically);
                 if non-NULL, user-given y-ticks (on original scale);
                 can be a list with one (numeric or NULL) item per panel}
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
  \item{withSubst}{logical; if \code{TRUE} (default) pattern substitution for
      titles and lables is used; otherwise no substitution is used. }
  \item{cex.pts}{size of the points of the second argument plotted, can be a vector;
         if it has the same length as the number of observations prior to any
         selection then it overrides argument \code{cex.nonlbl} and assigns 
         individual sizes to each of the points; otherwise it is recycled to
         the length of selected observations.}
  \item{col.pts}{color of the points of the second argument plotted, can
         be a vector as in \code{cex.pts}.}
  \item{pch.pts}{symbol of the points of the second argument plotted, can
         be a vector as in \code{cex.pts}.}
  \item{with.lab}{logical; shall labels be plotted to the observations?}
  \item{lab.pts}{character or NULL; labels to be plotted to the observations; if \code{NULL}
                 observation indices;}
  \item{lab.font}{font to be used for labels (of the observations).}       
  \item{alpha.trsp}{alpha transparency to be added ex post to colors
        \code{col.pch} and \code{col.lbl}; if one-dim and NA all colors are
        left unchanged. Otherwise, with usual recycling rules \code{alpha.trsp}
        gets shorted/prolongated to length the data-symbols to be plotted.
        Coordinates of this vector \code{alpha.trsp} with NA are left unchanged,
        while for the remaining ones, the alpha channel in rgb space is set
        to the respective coordinate value of \code{alpha.trsp}. The non-NA
        entries must be integers in [0,255] (0 invisible, 255 opaque).}
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
   otherwise we return \code{invisible()} as usual.}
  \item{\dots}{further parameters for \code{plot}}
}
\value{
  The function returns invisibly a list of elements containing the information 
  needed to compute the respective diagnostic plot.
   The return value allows to recover all information used to produce the plot
   for later use in enhanced graphics (e.g. with ggplot).   
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
has correct length. If argument \code{withSubst} is \code{TRUE}, in all title 
and axis lable arguments, the following patterns are substituted:
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

In addition, argument \code{\dots} may contain arguments \code{panel.first},
\code{panel.last}, i.e., hook expressions to be evaluated at the very beginning
and at the very end of each panel (within the then valid coordinates).
To be able to use these hooks for each panel individually, they may also be
lists of expressions (of the same length as the number of panels and
run through in the same order as the panels).
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
     ylim=c(-3,3,-1,3), xlim=c(-2,3),
     with.legend = TRUE)

data <- r(N)(30)
plot(IC2, data, panel.first= grid(),
     ylim = c(-3,3,-1,3), xlim=c(-2,3),
     cex.pts = 3, pch.pts = 1:2, col.pts="green",
     with.lab = TRUE, which.lbs = c(1:4,15:20),
     which.Order = 1:6, return.Order = TRUE)
}
\keyword{methods}
\keyword{distribution}

