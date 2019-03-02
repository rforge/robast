\name{ComparePlot}
\alias{ComparePlot}
\title{Wrapper function for function comparePlot}
\usage{
  ComparePlot(IC1, IC2,  y, ..., IC3 = NULL, IC4 = NULL,
    alpha.trsp = 100, with.legend = TRUE, rescale = FALSE,
    withCall = TRUE)
}
\arguments{
  \item{IC1}{object of class \code{IC}}

  \item{IC2}{object of class \code{IC}}

  \item{IC3}{object of class \code{IC}}

  \item{IC4}{object of class \code{IC}}

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
  \code{invisible(retV)} where \code{retV} is the return value
  of the respective call to the full-fledged function \code{comparePlot}
  with the additional item \code{wrapcall} with the call
  to the wrapper \code{ComparePlot} and \code{wrappedcall} the call to
  to the full-fledged function \code{comparePlot}.
}
\description{
  The wrapper \code{ComparePlot} (capital C!) takes most of arguments 
  to function \code{comparePlot} (lower case c!) by default and gives a user 
  possibility to run the function with low number of arguments.
}
\section{Details}{
  Calls \code{comparePlot} with suitably chosen defaults;
  if \code{withCall == TRUE}, the call to
  \code{comparePlot}, i.e., item \code{wrappedcall} of the (hidden)
  return value, is printed.
}
\examples{
# Gamma
fam <- GammaFamily()
rfam <- InfRobModel(fam, ContNeighborhood(0.5))
IC1 <- optIC(model = fam, risk = asCov())
IC2 <- makeIC(list(function(x)sin(x),function(x)x^2), L2Fam = fam)
Y <- distribution(fam)
y <- r(Y)(100)
ComparePlot(IC1, IC2, y, withCall = TRUE)
}

