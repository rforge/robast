\name{KolmogorovMinDist}
\alias{KolmogorovMinDist}
\alias{KolmogorovMinDist-methods}
\alias{KolmogorovMinDist,matrix,Norm-method}
\alias{KolmogorovMinDist,AffyBatch,AbscontDistribution-method}
\alias{KolmogorovMinDist,BeadLevelList,AbscontDistribution-method}

\title{Generic Function for Computing Minimum Kolmogorov Distance for Biological Data}
\description{
  Generic function for computing minimum Kolmogorov distance for biological data.
}
\usage{
KolmogorovMinDist(x, D, ...)

\S4method{KolmogorovMinDist}{matrix,Norm}(x, D, mad0 = 1e-4)

\S4method{KolmogorovMinDist}{AffyBatch,AbscontDistribution}(x, D, bg.correct = TRUE, pmcorrect = TRUE, 
                  verbose = TRUE)

\S4method{KolmogorovMinDist}{BeadLevelList,AbscontDistribution}(x, D, log = FALSE, imagesPerArray = 1, what = "G", 
                  probes = NULL, arrays = NULL)
}
\arguments{
  \item{x}{ biological data. }
  \item{D}{ object of class \code{"UnivariateDistribution"}. }
  \item{\dots}{ additional parameters. }
  \item{mad0}{ scale estimate used if computed MAD is equal to zero. 
    Median and MAD are used as start parameter for optimization. }
  \item{bg.correct}{ if \code{TRUE} MAS 5.0 background correction is performed;
    confer \code{\link[affy:bgc]{bg.correct.mas}}. }
  \item{pmcorrect}{ if \code{TRUE} \code{log2(PM/MM)} is used. 
    If \code{FALSE} only \code{log2(PM)} is used. }
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed. }
  \item{log}{ if \code{TRUE}, then the log2 intensities for each bead-type are summarized. }
  \item{imagesPerArray}{ Specifies how many images (strips) there are per array. 
    Normally 1 for a SAM and 1 or 2 for a BeadChip. The images (strips) from the same array 
    will be combined so that each column in the output represents a sample. }
  \item{what}{ character string specifying which intensities/values to summarize. 
    See \code{\link[beadarray]{getArrayData}} for a list of possibilities. }
  \item{probes}{ Specify particular probes to summarize. If left \code{NULL} then all
    the probes on the first array are used. }
  \item{arrays}{ integer (scalar or vector) specifying the strips/arrays to summarize. 
    If \code{NULL}, then all strips/arrays are summarized. }
}
\details{
  The minimum Kolmogorov distance is computed for each row of a matrix, each
  Affymetrix probe, or each Illumina bead, respectively.
  
  So far, only the minimum distance to the set of normal distributions can be
  computed.
}
\value{ List with components \code{dist} containing a numeric vector 
or matrix with minimum Kolmogorov distances and \code{n} a numeric vector 
or matrix with the corresponding sample sizes. }
\references{
  Huber, P.J. (1981) \emph{Robust Statistics}. New York: Wiley.

  Kohl M. and Deigner H.P. (2009). Using infinitesimally robust estimators for 
  preprocessing gene expression data. In preparation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[distrEx]{KolmogorovDist}}, \code{\link[distrMod]{MDEstimator}} }
\examples{
## matrix method for KolmogorovMinDist
ind <- rbinom(200, size=1, prob=0.05) 
X <- matrix(rnorm(200, mean=ind*3, sd=(1-ind) + ind*9), nrow = 2)
KolmogorovMinDist(X, D = Norm())

## using Affymetrix-Data
data(SpikeIn)
probes <- log2(pm(SpikeIn))
(res <- KolmogorovMinDist(probes, Norm()))
boxplot(res$dist)

\dontrun{
## "Not run" just because of computation time
require(affydata)
data(Dilution)
res <- KolmogorovMinDist(Dilution[,1], Norm())
summary(res$dist)
boxplot(res$dist)
plot(res$n, res$dist, pch = 20, main = "Kolmogorov distance vs. sample size",
     xlab = "sample size", ylab = "Kolmogorov distance",
     ylim = c(0, max(res$dist)))
uni.n <- min(res$n):max(res$n)
lines(uni.n, 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
}

## using Illumina-Data
\dontrun{
## "Not run" just because of computation time
data(BLData)
res <- KolmogorovMinDist(BLData, Norm(), arrays = 1)
res1 <- KolmogorovMinDist(BLData, log = TRUE, Norm(), arrays = 1)
summary(cbind(res$dist, res1$dist))
boxplot(list(res$dist, res1$dist), names = c("raw", "log-raw"))
sort(unique(res1$n))
plot(res1$n, res1$dist, pch = 20, main = "Kolmogorov distance vs. sample size",
     xlab = "sample size", ylab = "Kolmogorov distance",
     ylim = c(0, max(res1$dist)), xlim = c(min(res1$n), 56))
uni.n <- min(res1$n):56
lines(uni.n, 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
}
}
\concept{normal location and scale}
\concept{minimum distance}
\keyword{robust}
