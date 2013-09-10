\name{internals_for_RobAStBase_plot}
\alias{internals_for_RobAStBase_plot}
\alias{.rescalefct}
\alias{.plotRescaledAxis}
\alias{.legendCoord}
\alias{.SelectOrderData}
\alias{.makedotsP}
\alias{.makedotsLowLevel}

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
                  x.ticks = NULL, y.ticks = NULL, withbox = TRUE)
.legendCoord(x, scaleX, scaleX.fct, scaleY, scaleY.fct)
.SelectOrderData(data, fct, which.lbs, which.Order)
.makedotsP(dots)
.makedotsLowLevel(dots)
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
  \item{withbox}{logical of length 1. If \code{TRUE}, even if \code{scaleX} and
   \code{scaleY} are both \code{FALSE} and, simultaneously, \code{x.ticks} and
   \code{y.ticks} are both \code{NULL}, a respective box is drawn around the
   panel; otherwise no box is drawn in this case. }
  \item{dots}{a list; intended to be the \code{\dots} argument of \code{plot}. }
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
indices of the data selected by \code{which.lbs} in the original data.

\code{.makedotsP} and \code{.makedotsLowLevel} manipulate the \code{\dots}
argument, deleting certain items and selecting items which can be digested by
\code{plot}, returning the manipulated list.
}


\keyword{internal}
\concept{utilities}
\keyword{hplot}
