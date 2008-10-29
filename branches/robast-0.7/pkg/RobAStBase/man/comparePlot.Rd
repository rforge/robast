\name{comparePlot-methods}
\docType{methods}
\alias{comparePlot}
\alias{comparePlot-methods}
\alias{comparePlot,IC,IC-method}

\title{Compare - Plots}
\description{
  Plots 2-4 influence curves to the same model.
}
\usage{
comparePlot(obj1, obj2, ... )
\S4method{comparePlot}{IC,IC}(obj1, obj2, obj3 = NULL, obj4 = NULL, 
             ..., withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3],
             mfColRow = TRUE)
}
\arguments{
  \item{obj1}{ object of class \code{"InfluenceCurve"} }
  \item{obj2}{ object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{obj3}{ optional: object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{obj4}{ optional: object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
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
  \item{cex.inner}{magnification to be used for inner titles relative
          to the current setting of \code{cex}; as in 
          \code{\link[stats]{par}}}
  \item{col.inner}{character or integer code; color for the inner title}              
  \item{mfColRow}{shall default partition in panels be used --- defaults to \code{TRUE}}
  \item{\dots}{further arguments to be passed to \code{plot}}
}
\details{
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
\item{\code{"\%C1"},\code{"\%C2"},[\code{"\%C3"},][\code{"\%C4"}]}{class of argument 
      \code{obj<i>}, i=1,..4}
\item{\code{"\%A1"},\code{"\%A2"},[\code{"\%A3"},][\code{"\%A4"}]}{deparsed argument  
      \code{obj<i>}, i=1,..4}
\item{\code{"\%D"}}{time/date-string when the plot was generated}
}
}

%\value{}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}, \code{\link[base]{plot}}}
\examples{
if(require(ROptEst)){

N0 <- NormLocationScaleFamily(mean=0, sd=1) 
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))

IC1 <- optIC(model = N0, risk = asCov())
IC2 <- optIC(model = N0.Rob1, risk = asMSE())

comparePlot(IC1,IC2)
}
}
\keyword{robust}
