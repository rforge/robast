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
\S4method{infoPlot}{IC}(object, ..., withSweave = getdistrOption("withSweave"), 
             colI = grey(0.5), lwdI = 0.7*par("lwd"),
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3],
             mfColRow = TRUE)
}
\arguments{
  \item{object}{ object of class \code{"InfluenceCurve"} }
  \item{withSweave}{logical: if \code{TRUE} (for working with \command{Sweave}) 
        no extra device is opened}
  \item{main}{logical: is a main title to be used? or \cr
              just as argument \code{main} in \code{\link{plot.default}}.}
  \item{inner}{logical: panels have their own titles? or \cr
               character vector of / cast to length number of comparands: 
               \code{main} in \code{\link{plot.default}})} 
  \item{sub}{logical: is a sub-title to be used? or \cr
              just as argument \code{sub} in \code{\link{plot.default}}.}
  \item{tmar}{top margin -- useful for non-standard main title sizes}
  \item{bmar}{bottom margin -- useful for non-standard sub title sizes}
  \item{colI}{color of the classically optimal IC}
  \item{lwdI}{linewidth of the classically optimal IC}
  \item{cex.inner}{magnification to be used for inner titles relative
          to the current setting of \code{cex}; as in 
          \code{\link[stats]{par}}}
  \item{col.inner}{character or integer code; color for the inner title}              
  \item{mfColRow}{shall default partition in panels be used --- defaults to \code{TRUE}}
  \item{\dots} {further parameters for \code{plot}}
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
has correct length. In all title arguments, the following patterns are substituted:
\itemize{
\item{\code{"\%C"}}{class of argument \code{object}}
\item{\code{"\%A"}}{deparsed argument  \code{object}}
\item{\code{"\%D"}}{time/date-string when the plot was generated}
}
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
}
\concept{absolute information}
\concept{relative information}
\keyword{robust}
