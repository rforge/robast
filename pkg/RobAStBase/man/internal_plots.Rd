\name{internals_for_RobAStBase_plot}
\alias{internals_for_RobAStBase_plot}
\alias{.rescalefct}
\alias{.plotRescaledAxis}
\alias{.legendCoord}
\alias{.SelectOrderData}
\alias{.makedotsP}
\alias{.makedotsLowLevel}
\alias{.cexscale}
\alias{.getX.vec}
\alias{.getXlimYlim}
\alias{.prepareTitles}
\alias{.getToDraw}
\alias{.preparePanelFirstLast}

\title{Internal / Helper functions of package RobAStBase for plot functions}

\description{
These functions are internally used helper functions for \code{\link{plot}},
\code{\link{infoPlot}} \code{\link{comparePlot}} in package \pkg{RobAStBase}.}

\usage{
.rescalefct(x, fct, scaleX = FALSE, scaleX.fct, scaleX.inv,
         scaleY = FALSE, scaleY.fct = pnorm,
         xlim, ylim, dots)
.plotRescaledAxis(scaleX, scaleX.fct, scaleX.inv, scaleY,scaleY.fct,
                  scaleY.inv, xlim, ylim, X, ypts = 400, n = 11,
                  finiteEndpoints = rep(FALSE,4),
                  x.ticks = NULL, y.ticks = NULL, withbox = TRUE)
.legendCoord(x, scaleX, scaleX.fct, scaleY, scaleY.fct)
.SelectOrderData(data, fct, which.lbs, which.Order, which.nonlbs = NULL)
.makedotsP(dots)
.makedotsLowLevel(dots)
.cexscale(y, y1=y, maxcex=4,mincex=0.05,cex, fun=NULL)
.getX.vec(distr, dims0, lty, x.vec, scaleX, scaleX.fct, scaleX.inv, xm, xM)
.getXlimYlim(dots,dotsP, dims0, xlim, ylim)
.prepareTitles(withSubst, presubArg2, presubArg3, dots, mainText,
               L2Fam, inner, dims0, dims, to.draw, trafO, obj, type, bmar, tmar)
.getToDraw(dims, trafO, L2Fam, to.draw.arg, Abs=NULL)
.preparePanelFirstLast(with.automatic.grid , dims0, pF.0, pL.0,
            logArg, scaleX, scaleY, x.ticks, y.ticks, scaleX.fct, scaleY.fct)
}
\arguments{
  \item{scaleX}{logical; shall X-axis be rescaled (by default according to the cdf of
          the underlying distribution)?}
  \item{scaleY}{logical; shall Y-axis be rescaled (by default according to a probit scale)?}
  \item{scaleX.fct}{an isotone, vectorized function mapping the domain of the IC
            to [0,1]; if \code{scaleX} is \code{TRUE} and \code{scaleX.fct} is
            missing, the cdf of the underlying observation distribution.}
  \item{scaleX.inv}{the inverse function to \code{scale.fct}, i.e., an isotone,
            vectorized function mapping [0,1] to the domain of the IC
            such that for any \code{x} in the domain,\cr
            \code{scaleX.inv(scaleX.fct(x))==x}; if \code{scaleX} is \code{TRUE}
            and \code{scaleX.inv} is
            missing, the quantile function of the underlying observation distribution.}
  \item{scaleY.fct}{an isotone, vectorized function mapping for each coordinate the
            range of the respective coordinate of the IC
            to [0,1]; defaulting to the cdf of \eqn{{\cal N}(0,1)}{N(0,1)}.}
  \item{scaleY.inv}{an isotone, vectorized function mapping for each coordinate
            the range [0,1] into the range of the respective coordinate of the IC;
            defaulting to the quantile function of  \eqn{{\cal N}(0,1)}{N(0,1)}.}
  \item{fct}{for \code{.rescalefct}: (vectorized) function from reals to reals: function to produce y values.
  for \code{.SelectOrderData}: (vectorized) function from data to reals:
    function to produce the ordering.}
  \item{xlim}{numeric vector of length 2: limits of the plotted x region (in original
              scale).}
  \item{ylim}{numeric vector of length 2: limits of the plotted y region (in original
              scale).}
  \item{finiteEndpoints}{a logical of length 4: are the  unscaled \code{xlim[1]},
                         \code{xlim[2]}, \code{ylim[1]}, \code{ylim[2]} finite? }
  \item{x.ticks}{numeric: coordinates in original scale of user-given ticks on x-axis.}
  \item{y.ticks}{numeric: coordinates in original scale of user-given ticks on y-axis.}
  \item{n}{integer: number of default ticks in x and y axis.}
  \item{X}{numeric: x-values on [0,1]-scale from which to produce range of x axis.}
  \item{ypts}{integer: number of y-points to be generated on [0,1] range to
              produce y-axis.}
  \item{dots}{dots argument for later calls to \code{plot} which is manipulated
           when using rescaling.}
  \item{x}{coordinates of the legend placement to be rescaled when using a
      scale transformation.}
  \item{data}{numeric or data.frame containing numerics; observations to be plotted into
     a graph.}
  \item{which.lbs}{either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then no observation is excluded}
  \item{which.Order}{indices of the observations after ordering which are to be
    plotted into the graph; more speficically, the indices are obtained after
    a possible thin-out by \code{which.lbs} and after ordering in descending order
    of the remaining observations. If this argument is \code{NULL} then no
    (further) observation is excluded.}
  \item{which.nonlbs}{indices of the observations which should be plotted but
        not labelled; either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then all non-labelled
          observations are plotted.}
  \item{withbox}{logical of length 1. If \code{TRUE}, even if \code{scaleX} and
   \code{scaleY} are both \code{FALSE} and, simultaneously, \code{x.ticks} and
   \code{y.ticks} are both \code{NULL}, a respective box is drawn around the
   panel; otherwise no box is drawn in this case. }
  \item{dots}{a list; intended to be the \code{\dots} argument of \code{plot}. }
  \item{y}{input at which to compute \code{cex}-values (for these \code{y} values)}
  \item{y1}{secondary input of \code{y}-values which is also used for
         computaiton but not plotted itself at this \code{cex} value.}
  \item{maxcex}{maximal value for \code{cex}}
  \item{mincex}{minimal value for \code{cex}}
  \item{cex}{basic factor for \code{cex}}
  \item{fun}{(optional) function to compute the \code{cex}-values acc. to the
              \code{y} and \code{y1}-values.}
\item{L2Fam}{L2-Family},
\item{to.draw.arg}{resp. argument from the user-interface of \code{comparePlot}, \code{infoPlot}, and the
\code{plot}-method for ICs, i.e., either \code{NULL} (default;
everything is plotted) or a vector of either integers (the indices of the
subplots to be drawn) or characters --- the names of the subplots to be drawn:
these names are to be chosen either among the row names of the trafo matrix
\code{rownames(trafo(eval(x@CallL2Fam)@param))} or if the last expression is
\code{NULL} a vector \code{"dim<dimnr>"}, \code{dimnr} running through
the number of rows of the trafo matrix.}
\item{panelFirst}{the \code{panel.first} expression to be extended}
\item{distr}{underlying distribution of the observations}
\item{dims0}{number of panels to be plotted}
\item{lty}{the \code{lty} argument used in the user interface}
\item{x.vec}{the resp. argument from the user interface.}
\item{xm}{minimal x-value}
\item{xM}{maximal x-value}
\item{dots}{the internal \code{dots}-variable to be changed}
\item{dotsP}{the internal \code{dotsP}-variable to be changed (taking up
  additional arguments for the calls to \code{plot})}
\item{xlim}{the possibly matrix-valued argument of the user interface}
\item{ylim}{the possibly matrix-valued argument of the user interface}
\item{withSubst}{should substitutions be made in titles?}
\item{presubArg2}{title part}
\item{presubArg3}{title part}
\item{mainText}{text to be used for \code{main}}
\item{inner}{the resp. argument from the user interface.}
\item{dims}{the dimension of the pIC}
\item{to.draw}{the actual numbers / names of the plotted panels.}
\item{trafO}{the parameter transformation of  \code{L2fam@param}.}
\item{obj}{the object from which to draw the name (of the plot).}
\item{type}{for which type of plot should the titles be produced -- \code{"info"},
\code{"all"} or \code{compare}.}
\item{bmar}{resp. argument from the user-interface of \code{comparePlot}, \code{infoPlot}, and the
\code{plot}-method for ICs.}
\item{tmar}{resp. argument from the user-interface of \code{comparePlot}, \code{infoPlot}, and the
\code{plot}-method for ICs.}
\item{Abs}{additional text template used in \code{infoPlot}.}
\item{with.automatic.grid}{}
\item{pF.0}{the \code{panel.first} expression given by rhe user}
\item{pL.0}{the \code{panel.last} expression  given by rhe user}
\item{logArg}{a character vector with \code{log}-arguments for  \code{plot}
for each of the plotted panels.}
}
\details{
\code{.rescalefct} rescales, if necessary, x and y axis for use in plot
   functions. More specifically, if \code{scaleX} is \code{TRUE} rescales x,
   if \code{scaleY} is \code{TRUE} rescales fct(x) (otherwise leaves them
   unchanges); to this end uses trafos \code{scaleX.fct} with inverse
   \code{scaleX.inv}, resp. \code{scaleY.fct}; it respects \code{xlim}
   and  \code{ylim} (given in orig. scale), thins out the scaled values if
   necessary and accordingly modifies  slots \code{xaxt}, \code{yaxt}, and
   \code{axes} of argument \code{dots} to indicate the new axes have to be drawn;
   using the  paradigm small letters to denote values on original scale and
    capital letters on transformed scale, its return value is a list with
    (thinned out) values of x and y, X and Y and modified dots.

\code{.plotRescaledAxis} plots rescaled axes according to logicals
\code{scaleX}, \code{scaleY}; to this end uses trafos \code{scaleX.fct} with
inverse \code{scale.inv}, resp. \code{scaleY.fct}, \code{scaleY.inv};
 it respects \code{xlim} and \code{ylim}. By default, ot produces the x axes according
 to the values in argument \code{X}, and the y axes as an equidistant grid of
 length \code{ypts} on [0,1] (on transformed scale); each of these axes, again
 by default will have \code{n} tick values; these are however thinned out if
 the come to lie too close to each other on transformed scale.
 Instead of producing automatically chosen tick values, the user may explicitly
 require x-ticks and y-ticks values on the axes, using arguments \code{x.ticks}
 and \code{y-ticks}. This function has no return value.

\code{.legendCoord} produces, if needed (i.e., if coordinates are not
given as strings like \code{"bottomright"}), rescaled coordinates for the
placement of a legend.

\code{.SelectOrderData}, for data to be plotted into the graph,
 performs two optional selections: a first selection on the unordered (original)
 data (acc. to argument \code{which.lbs}) and a second selection according to
 \code{which.Order} on the data remaining after the first selection and ordered
 according to argument \code{fct}; the return value is a list with elements
\code{data}, ie., the selected/thinned out data, \code{y}, ie., the
values of  \code{fct(data)}, \code{ind}, ie., the indices of the selected data
in the original data (after possibly two selections), and \code{ind1} the
indices of the data selected by \code{which.lbs} in the original data; in 
addition also the non selected data, \code{data.ns}, the respective
y-values \code{y.ns} and the corresponding index elements \code{ind.ns} are
returned as list items.

\code{.makedotsP} and \code{.makedotsLowLevel} manipulate the \code{\dots}
argument, deleting certain items and selecting items which can be digested by
\code{plot}, returning the manipulated list.

\code{.cexscale}{rescales the point sizes of the points to be plotted;
     the unscaled sizes are given in argument \code{y}, \code{y1} in
     case of several lines of points to be plotted may contain the
     vector of the sizes of all points to be plotted in (e.g., including
     those of the other lines of points). \code{maxcex} and \code{mincex}
     are maximum and minimum of the raw rescaled sizes; \code{cex} is
     a factor drawn from argument \code{cex.pts} by which the raw sizes
     are rescaled before being returned. \code{fun} is the function
     by which the rescaling is done; by default this argument is \code{NULL}
     and in this case the function \code{log(1+abs(x))} is used.
}
\code{.getDimsTD}{returns the number of different coordinates to be plotted.}
\code{.producePanelFirstS}{for each graphical panel inserts
                           (if needed) x and y tickmarks for user-specific axes
                           into a \code{panel.first} expression.}
\code{.getX.vec}{produces the x-grid (on original scale) for each of the panels
                 and decides whether to plot lines or points.}

\code{.getXlimYlim}{produces panel-wise \code{xlim} and \code{ylim} arguments
                    as well as left and right endpoints of the x-scalas.}
\code{.prepareTitles}{produces the titles for the panels}
\code{.getToDraw}{computes which panels are to be drawn}
\code{.preparePanelFirstLast}{prepares the \code{panel.first} and \code{panel.last}
      expressions for each of the panels.}
}


\keyword{internal}
\concept{utilities}
\keyword{hplot}
