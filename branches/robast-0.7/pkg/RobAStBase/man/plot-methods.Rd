\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,IC,missing-method}
\title{ Methods for Function plot in Package `RobAStBase' }
\description{plot-methods}
\usage{plot(x, y, ...)
\S4method{plot}{IC,missing}(x, ..., withSweave = getdistrOption("withSweave"),
             main = FALSE, inner = TRUE, sub = FALSE,
             col.inner = par("col.main"), cex.inner = 0.8,
             bmar = par("mar")[1], tmar = par("mar")[3],
             mfColRow = TRUE)
}
\arguments{
  \item{x}{object of class \code{"IC"}: IC to be plotted }
  \item{y}{ missing }
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
  \item{\dots} {further parameters for \code{plot}}
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
\item{\code{"\%C"}}{class of argument \code{object}}
\item{\code{"\%A"}}{deparsed argument  \code{object}}
\item{\code{"\%D"}}{time/date-string when the plot was generated}
}
}
\examples{
IC1 <- new("IC")
plot(IC1)
plot(IC1, main = TRUE, panel.first= grid(),
     col = "blue", cex.main = 2, cex.inner = 1)
}
\keyword{methods}
\keyword{distribution}

