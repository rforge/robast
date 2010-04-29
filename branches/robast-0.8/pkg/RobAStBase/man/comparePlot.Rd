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
             data = NULL, ..., withSweave = getdistrOption("withSweave"),
             main = FALSE, inner = TRUE, sub = FALSE, 
             col = par("col"), lwd = par("lwd"), lty, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3],
             with.legend = TRUE, legend.bg = "white",
             legend.location = "bottomright", legend.cex = 0.8,
             mfColRow = TRUE, to.draw.arg = NULL,
             cex.pts = 1, col.pts = par("col"),
             pch.pts = 1, jitter.fac = 1, with.lab = FALSE,
             lab.pts = NULL, lab.font = NULL,
             which.lbs = NULL, which.Order  = NULL, return.Order = FALSE)
}
\arguments{
  \item{obj1}{ object of class \code{"InfluenceCurve"} }
  \item{obj2}{ object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{obj3}{ optional: object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{obj4}{ optional: object of class \code{"InfluenceCurve"} to be compared with \code{obj1}}
  \item{data}{optional data argument --- for plotting observations into the plot;}
  \item{withSweave}{logical: if \code{TRUE} (for working with \command{Sweave})
        no extra device is opened}
  \item{main}{logical: is a main title to be used? or \cr
              just as argument \code{main} in \code{\link[graphics]{plot.default}}.}
  \item{col}{color[s] of ICs in arguments \code{obj1} [,\ldots,\code{obj4}].}
  \item{lwd}{linewidth[s] of ICs in arguments \code{obj1} [,\ldots,\code{obj4}].}
  \item{lty}{line-type[s] of ICs in arguments \code{obj1} [,\ldots,\code{obj4}].}
  \item{inner}{logical: do panels have their own titles? or \cr
               character vector of / cast to length 'number of plotted 
               dimensions';
               if argument \code{to.draw.arg} is used, this refers to 
               a vector of length \code{length(to.draw.arg)}, the 
               actually plotted dimensions. For further information, see also 
               description of argument \code{main} in \code{\link[graphics]{plot.default}}.} 
  \item{sub}{logical: is a sub-title to be used? or \cr
              just as argument \code{sub} in \code{\link[graphics]{plot.default}}.}
  \item{tmar}{top margin -- useful for non-standard main title sizes}
  \item{bmar}{bottom margin -- useful for non-standard sub title sizes}
  \item{cex.inner}{magnification to be used for inner titles relative
          to the current setting of \code{cex}; as in 
          \code{\link[graphics]{par}}}
  \item{col.inner}{character or integer code; color for the inner title}              
  \item{with.legend}{logical; shall a legend be plotted?}
  \item{legend.location}{a valid argument \code{x} for \code{\link{legend}} ---
                         the place where to put the legend on the last issued
                         plot}
  \item{legend.bg}{background color for the legend}
  \item{legend.cex}{magnification factor for the legend}
  \item{mfColRow}{shall default partition in panels be used --- defaults to \code{TRUE}}
  \item{to.draw.arg}{Either \code{NULL} (default; 
                          everything is plotted) or a vector of either integers 
                         (the indices of the subplots to be drawn) or characters 
                         --- the names of the subplots to be drawn: these
                         names are to be chosen either among the row names of 
                         the trafo matrix 
                         \code{rownames(trafo(eval(obj1@CallL2Fam)@param))} 
                         or if the last expression is \code{NULL} a 
                         vector \code{"dim<dimnr>"}, \code{dimnr} running through 
                         the number of rows of the trafo matrix.
                         }
  \item{cex.pts}{size of the points of the \code{data} argument plotted}
  \item{col.pts}{color of the points of the \code{data} argument plotted}
  \item{pch.pts}{symbol of the points of the \code{data} argument plotted}
  \item{with.lab}{logical; shall labels be plotted to the observations?}
  \item{lab.pts}{character or NULL; labels to be plotted to the observations; if \code{NULL}
                 observation indices;}
  \item{lab.font}{font to be used for labels}
  \item{jitter.fac}{jittering factor used in case of a \code{DiscreteDistribution}
                    for plotting points of the \code{data} argument in a jittered fashion.}
  \item{which.lbs}{either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then no observation is excluded}
  \item{which.Order}{for each of the given ICs, we order the observations (descending)
           according to the norm given by the corresponding \code{normtype(object)};
           then \code{which.Order} either is an integer vector with the indices of the \emph{ordered}
           observations (remaining after a possible reduction by argument \code{which.lbs})
           to be plotted into graph or \code{NULL} --- then no (further) observation
           is excluded.}
  \item{return.Order}{logical; if \code{TRUE}, a list of length maximally four
    with order vectors is returned --- one for the ordering w.r.t. each of
    the given ICs; more specifically, the order of the (remaining) observations
    given by their original index is returned (remaining means: after a possible
    reduction by argument \code{which.lbs}, and ordering is according to the norm given by
           \code{normtype(object)});
   othervise we return \code{invisible()} as usual.}
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
\describe{
\item{\code{"\%C1"},\code{"\%C2"},[\code{"\%C3"},][\code{"\%C4"}]}{class of argument 
      \code{obj<i>}, i=1,..4}
\item{\code{"\%A1"},\code{"\%A2"},[\code{"\%A3"},][\code{"\%A4"}]}{deparsed argument  
      \code{obj<i>}, i=1,..4}
\item{\code{"\%D"}}{time/date-string when the plot was generated}
}

If argument \code{\dots} contains argument \code{ylim}, this may either be
as in \code{plot.default} (i.e. a vector of length 2) or a vector of 
length 2*(number of plotted dimensions); in the case of longer length, 
these are the values for \code{ylim} for the plotted dimensions of the IC, 
one pair for each dimension.
}

%\value{}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.
}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link[distrMod]{L2ParamFamily-class}}, \code{\link{IC-class}}, \code{\link[graphics]{plot}}}
\examples{
if(require(ROptEst)){

N0 <- NormLocationScaleFamily(mean=0, sd=1) 
N0.Rob1 <- InfRobModel(center = N0, neighbor = ContNeighborhood(radius = 0.5))

IC1 <- optIC(model = N0, risk = asCov())
IC2 <- optIC(model = N0.Rob1, risk = asMSE())

comparePlot(IC1,IC2)

data <- r(N0)(20)
comparePlot(IC1, IC2, data=data, with.lab = TRUE,
            which.lbs = c(1:4,15:20),
            which.Order = 1:6,
            return.Order = TRUE)

## selection of subpanels for plotting
par(mfrow=c(1,1))
comparePlot(IC1, IC2 ,mfColRow = FALSE, to.draw.arg=c("mean"),
            panel.first= grid(),ylim=c(-4,4),xlim=c(-6,6))
## matrix-valued ylim
comparePlot(IC1, IC2, panel.first= grid(),ylim=c(-4,4,0,4),xlim=c(-6,6))

## with use of trafo-matrix:
G <- GammaFamily(scale = 1, shape = 2)
## explicitely transforming to
## MASS parametrization:
mtrafo <- function(x){
     nms0 <- names(c(main(param(G)),nuisance(param(G))))
     nms <- c("shape","rate")
     fval0 <- c(x[2], 1/x[1])
     names(fval0) <- nms
     mat0 <- matrix( c(0, -1/x[1]^2, 1, 0), nrow = 2, ncol = 2,
                     dimnames = list(nms,nms0))                          
     list(fval = fval0, mat = mat0)}
G2 <- G
trafo(G2) <- mtrafo
G2
G2.Rob1 <- InfRobModel(center = G2, neighbor = ContNeighborhood(radius = 0.5))
system.time(IC1 <- optIC(model = G2, risk = asCov()))
system.time(IC2 <- optIC(model = G2.Rob1, risk = asMSE()))
system.time(IC2.i <- optIC(model = G2.Rob1, risk = asMSE(normtype=InfoNorm())))
system.time(IC2.s <- optIC(model = G2.Rob1, risk = asMSE(normtype=SelfNorm())))

comparePlot(IC1,IC2, IC2.i, IC2.s)


}
}
\keyword{robust}
