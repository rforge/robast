\name{internals_for_ROptEst_Cniperplots}
\alias{internals_for_ROptEst_Cniperplots}
\alias{.plotData}
\alias{.getFunCnip}
\title{Internal / Helper functions of package ROptEst for cniper plot functions}

\description{
These functions are internally used helper function for \code{\link{cniperCont}}
and \code{cniperPointPlot}.}

\usage{
.plotData(data, dots, fun, L2Fam,  IC, jit.fac, jit.tol, plotInfo )
.getFunCnip(IC1,IC2, risk, L2Fam, r, b20=NULL)
}
\arguments{
  \item{data}{data to be plot in}
  \item{dots}{list; argument \code{\dots} from the calling function}
  \item{origCl}{call: call from the calling function.}
  \item{fun}{function from data to reals; function according to which
  the data is ordered}
  \item{L2Fam}{ object of class \code{L2ParamFamily} }
  \item{IC}{object of class \code{IC}: To determine the size of the observation
    symbols, we compute the norm of the IC evaluated at this point.}
  \item{IC1}{object of class \code{IC}: one of the ICs from which to take the
  risk difference.}
  \item{IC2}{object of class \code{IC}: one of the ICs from which to take the
  risk difference.}
  \item{risk}{object of class \code{asGRisk}: the risk to take the fifference
  from}
  \item{r}{positive numeric of length 1: the neighborhood radius}
  \item{b20}{positive numeric of length 1: the maximal bias of \code{IC2}.}
  \item{jit.fac}{jittering factor used in case of a \code{DiscreteDistribution}
                    for plotting points of the second argument in a jittered fashion.}
  \item{jit.tol}{jittering tolerance used in case of a \code{DiscreteDistribution}
                    for plotting points of the second argument in a jittered fashion.}
  \item{plotInfo}{stored info from the plot}
}
\details{
\code{.plotData} takes argument \code{data} and plots it into the cniper graph.

\code{.getFunCnip} produces a function to compute the risk difference. If
argument \code{b20} is not \code{NULL}, in the risk difference, for \code{IC2}
     uses the least favorable contamination situation
     ('over all real Dirac contamination points'), i.e. leading to a bias
     of \code{b20}. Otherwise it uses the bias obtaine from a contamination
     in the evaluation point.

}
\value{
\item{.plotData}{
If argument \code{dots} contains an element \code{return.Order} and
if this equals \code{TRUE}, the order of the (remaining) observations
given by their original index is returned (remaining means: after a possible
reduction by argument \code{which.lbs}, and ordering is according \code{fun};
otherwise we return \code{invisible()} as usual.}
\item{.getFunCnip}{
  a vectorized function to compute the risk difference.
}
}

\keyword{internal}
\concept{utilities}
\keyword{hplot}
