\name{internal_helpers_for_producing_grids_in_plots_RobAStBase}
\alias{internal_helpers_for_producing_grids_in_plots_RobAStBase}
\alias{.getDimsTD}
\alias{.producePanelFirstS}
\alias{.producePanelFirstSn}

\title{Internal / Helper functions of package RobAStBase for grids in plot functions}

\description{
These functions are internally used helper functions for \code{\link{plot}},
\code{\link{infoPlot}} \code{\link{comparePlot}} in package \pkg{RobAStBase}.}

\usage{
.getDimsTD(L2Fam,to.draw.arg)
.producePanelFirstS(panelFirst,IC,to.draw.arg, isInfoPlot=FALSE,
                                x.ticks, scaleX, scaleX.fct,
                                y.ticks, scaleY, scaleY.fct)
.producePanelFirstSn(panelFirst, x.ticks, scaleX, scaleX.fct,
                                y.ticks, scaleY, scaleY.fct, logArg)
}
\arguments{
  \item{L2Fam}{the model at which the plot is produced (of class \code{L2ParamFamily}).}
  \item{panelFirst}{argument \code{panel.first} to be mingled for grid plotting.}
  \item{x.ticks}{numeric: coordinates in original scale of user-given ticks on x-axis.}
  \item{scaleX}{logical; shall X-axis be rescaled (by default according to the cdf of
          the underlying distribution)?}
  \item{scaleX.fct}{an isotone, vectorized function mapping the domain of the IC
            to [0,1]; if \code{scaleX} is \code{TRUE} and \code{scaleX.fct} is
            missing, the cdf of the underlying observation distribution.}
  \item{y.ticks}{numeric: coordinates in original scale of user-given ticks on y-axis.}
  \item{scaleY}{logical; shall Y-axis be rescaled (by default according to a probit scale)?}
  \item{scaleY.fct}{an isotone, vectorized function mapping for each coordinate the
            range of the respective coordinate of the IC
            to [0,1]; defaulting to the cdf of \eqn{{\cal N}(0,1)}{N(0,1)}.}
  \item{logArg}{a character string which contains "x" if the x axis is to be logarithmic,
            "y" if the y axis is to be logarithmic and "xy" or "yx" if both axes are
            to be logarithmic; see \code{\link{plot.default}}.}
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
  \item{IC}{object of class \code{"InfluenceCurve"} }
  \item{isInfoPlot}{logical; is this function to be used in \code{infoPlot} or
                    (\code{TRUE}) in another plot (\code{FALSE})? }
}
\details{
\code{.getDimsTD} computes the number of panels to be plotted.
\code{.producePanelFirstS} produces an unevaluated expression to be
used as argument \code{panel.first} in the diagnostic plots; i.e.;
knowing the actual tickmarks of the axis at the time of evaluation,
code is inserted to plot horizontal and vertical grid lines through
these tickmarks.
}


\keyword{internal}
\concept{utilities}
\keyword{hplot}
