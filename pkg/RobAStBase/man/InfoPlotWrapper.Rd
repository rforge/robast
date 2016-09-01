\name{InfoPlot}
\alias{InfoPlot}
\title{Wrapper function for information plot method}
\usage{
  InfoPlot(IC, data, ..., alpha.trsp = 100,
    with.legend = TRUE, rescale = FALSE, withCall = TRUE)
}
\arguments{
  \item{IC}{object of class \code{IC}}

  \item{data}{optional data argument --- for plotting
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
  The wrapper \code{InfoPlot} (captial I!) takes most of arguments to the plot 
  method \code{infoPlot} (lower case i!) by default and gives a user possibility to run 
  the function with low number of arguments.
}
\section{Details}{
  Calls \code{infoPlot} with suitably chosen defaults. If
  \code{withCall == TRUE}, the call to \code{infoPlot} is
  returned
}
\examples{
# Gamma
fam  <-  GammaFamily()
IC <- optIC(model = fam, risk = asCov())
Y <- distribution(fam)
data  <-  r(Y)(500)
InfoPlot(IC, data, withCall = FALSE)
}

