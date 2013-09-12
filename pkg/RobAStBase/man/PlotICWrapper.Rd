\name{PlotIC}
\alias{PlotIC}
\title{Wrapper function for plot method for IC}
\usage{
  PlotIC(IC, y, ..., alpha.trsp = 100, with.legend = TRUE,
    rescale = FALSE, withCall = TRUE)
}
\arguments{
  \item{IC}{object of class \code{IC}}

  \item{y}{optional data argument --- for plotting
  observations into the plot}

  \item{...}{additional parameters (in particular to be
  passed on to \code{plot})}

  \item{alpha.trsp}{the transparency argument (0 to 100)
  for ploting the data}

  \item{with.legend}{the flag for showing the legend of the
  plot}

  \item{rescale}{the flag for rescaling the axes for better
  view of the plot}

  \item{withCall}{the flag for the call output}
}
\value{
  invisible(NULL)
}
\description{
  The wrapper \code{PlotIC} takes most of arguments to the \code{plot} 
  method by default and gives a user possibility to run the function
  with low number of arguments.
}
\section{Details}{
  Calls \code{plot} with suitably chosen defaults; if
  \code{withCall == TRUE}, the call to \code{plot} is
  returned.
}
\examples{
# Gamma
fam <- GammaFamily()
rfam <- InfRobModel(fam, ContNeighborhood(0.5))
IC <- optIC(model = fam, risk = asCov())
Y <- distribution(fam)
y <- r(Y)(1000)
PlotIC(IC, y, withCall = FALSE)
}

