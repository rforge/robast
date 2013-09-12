\name{CniperPointPlot}
\alias{CniperPointPlot}
\title{Wrapper function for cniperPointPlot - Computation and Plot
 of Cniper Contamination and Cniper Points}
\usage{
  CniperPointPlot(fam, ...,
    lower = getdistrOption("DistrResolution"),
    upper = 1 - getdistrOption("DistrResolution"),
    with.legend = TRUE, rescale = FALSE, withCall = TRUE)
}
\arguments{
  \item{fam}{object of class L2ParamFamily}

  \item{...}{additional parameters (in particular to be
  passed on to \code{plot})}

  \item{lower}{the lower end point of the contamination
  interval}

  \item{upper}{the upper end point of the contamination
  interval}

  \item{with.legend}{the flag for showing the legend of the
  plot}

  \item{rescale}{the flag for rescaling the axes for better view of the plot}

  \item{withCall}{the flag for the call output}
}
\value{
  invisible(NULL)
}
\description{
  The wrapper \code{CniperPointPlot} (capital C!) takes most of arguments 
  to the \code{cniperPointPlot} (lower case c!) function by default and gives 
  a user possibility to run the function with low number of arguments.
}
\section{Details}{
  Calls \code{cniperPointPlot} with suitably chosen
  defaults; if \code{withCall == TRUE}, the call to
  \code{cniperPointPlot} is returned.
}
\examples{
L2fam <- GammaFamily()
CniperPointPlot(fam=L2fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)
}

