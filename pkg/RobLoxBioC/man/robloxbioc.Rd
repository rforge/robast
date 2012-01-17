\name{robloxbioc}
\alias{robloxbioc}
\alias{robloxbioc-methods}
\alias{robloxbioc,matrix-method}
\alias{robloxbioc,AffyBatch-method}
\alias{robloxbioc,BeadLevelList-method}

\title{Generic Function for Preprocessing Biological Data}
\description{
  Generic function for preprocessing biological data using optimally robust
  (rmx) estimators; confer Rieder (1994), Kohl (2005), Rieder et al (2008).
}
\usage{
robloxbioc(x, ...)

\S4method{robloxbioc}{matrix}(x, eps = NULL, eps.lower = 0, eps.upper = 0.05, steps = 3L, 
           fsCor = TRUE, mad0 = 1e-4)

\S4method{robloxbioc}{AffyBatch}(x, bg.correct = TRUE, pmcorrect = TRUE, normalize = FALSE, 
           add.constant = 32, verbose = TRUE, eps = NULL, 
           eps.lower = 0, eps.upper = 0.05, steps = 3L, fsCor = TRUE, 
           mad0 = 1e-4, contrast.tau = 0.03, scale.tau = 10, 
           delta = 2^(-20), sc = 500)

\S4method{robloxbioc}{BeadLevelList}(x, log = TRUE, imagesPerArray = 1, what = "G", probes = NULL, 
           arrays = NULL, eps = NULL, eps.lower = 0, eps.upper = 0.05, 
           steps = 3L, fsCor = TRUE, mad0 = 1e-4)
}
\arguments{
  \item{x}{ biological data. }
  \item{\dots}{ additional parameters. }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
        See details below. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors. See details below. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors. See details below. }
  \item{steps}{ positive integer. k-step is used to compute the optimally robust estimator. }
  \item{fsCor}{ logical: perform finite-sample correction. See function \code{\link[RobLox]{finiteSampleCorrection}}. }
  \item{mad0}{ scale estimate used if computed MAD is equal to zero}
  \item{bg.correct}{ if \code{TRUE} MAS 5.0 background correction is performed;
    confer \code{\link[affy:bgc]{bg.correct.mas}}. }
  \item{pmcorrect}{ method used for PM correction; \code{TRUE} calls an algorithm which is 
    comparable to the algorithm of MAS 5.0; confer \code{\link[affy:pmcorrect]{pmcorrect.mas}}. 
    If \code{FALSE} only the PM intensities are used. }
  \item{normalize}{ logical: if \code{TRUE}, Affymetrix scale normalization is performed. }
  \item{add.constant}{ constant added to the MAS 5.0 expression values before the normalization
    step. Improves the variance of the measure one no longer devides by numbers close to 0
    when computing fold-changes. }
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed. }
  \item{contrast.tau}{ a number denoting the contrast tau parameter; confer the MAS 5.0 
    PM correction algorithm. }
  \item{scale.tau}{ a number denoting the scale tau parameter; confer the MAS 5.0 
    PM correction algorithm. }
  \item{delta}{ a number denoting the delta parameter; confer the MAS 5.0 
    PM correction algorithm. }
  \item{sc}{ value at which all arrays will be scaled to. }
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
  The optimally-robust resp. the radius-minimax (rmx) estimator for normal location 
  and scale is used to preprocess biological data. The computation uses a k-step 
  construction with median and MAD as starting estimators; cf. Rieder (1994) and 
  Kohl (2005).

  If the amount of gross errors (contamination) is known, it can be 
  specified by \code{eps}. The radius of the corresponding infinitesimal 
  contamination neighborhood (infinitesimal version of Tukey's gross error model) 
  is obtained by multiplying \code{eps} by the square root of the sample size. 

  If the amount of gross errors (contamination) is unknown, which is typically
  the case, try to find a rough estimate for the amount of gross errors, such that 
  it lies between \code{eps.lower} and \code{eps.upper}.

  If \code{eps} is \code{NULL}, the radius-minimax (rmx) estimator in sense of 
  Rieder et al. (2001, 2008), respectively Section 2.2 of Kohl (2005) is used.
  
  The algorithm used for Affymetrix data is similar to MAS 5.0 (cf. Affymetrix (2002)).
  The main difference is the substitution of the Tukey one-step estimator by our rmx 
  k-step (k >= 1) estimator in the PM/MM correction step. The optional scale normalization 
  is performed as given in Affymetrix (2002).
  
  In case of Illumina data, the rmx estimator is used to summarize the bead types. 
  The implementation for the most part was taken from function 
  \code{\link[beadarray]{createBeadSummaryData}}.

  For sample size <= 2, median and MAD are used for estimation.
  
  If \code{eps = 0}, mean and sd are computed.
}
\value{ Return value depends on the class of \code{x}. 
  In case of \code{"matrix"} a matrix with columns "mean" and "sd" is returned.
  In case of \code{"AffyBatch"} an object of class \code{"ExpressionSet"} is returned. 
}
\references{
  Affymetrix, Inc. (2002). \emph{Statistical Algorithms Description Document}.
  Affymetrix, Santa Clara.

  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \emph{17}(1) 13-40.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Submitted. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \url{www.uni-bayreuth.de/departments/math/org/mathe7/RIEDER/pubs/RR.pdf}
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[RobLox]{roblox}}, \code{\link[RobLox]{rowRoblox}},
         \code{\link[affy]{AffyBatch-class}}, 
         \code{\link[affy:generateExprVal.method.avgdiff]{generateExprVal.method.mas}},
         \code{\link[Biobase:class.ExpressionSet]{ExpressionSet-class}},
         \code{\link[beadarray]{createBeadSummaryData}} }
\examples{
set.seed(123) # to have reproducible results for package checking

## similar to rowRoblox of package RobLox
ind <- rbinom(200, size=1, prob=0.05)
X <- matrix(rnorm(200, mean=ind*3, sd=(1-ind) + ind*9), nrow = 2)
robloxbioc(X)
robloxbioc(X, steps = 5)
robloxbioc(X, eps = 0.05)
robloxbioc(X, eps = 0.05, steps = 5)

## the function is designed for large scale problems
X <- matrix(rnorm(50000*20, mean = 1), nrow = 50000)
system.time(robloxbioc(X))

## using Affymetrix-Data
## confer example to generateExprVal.method.mas
## A more worked out example can be found in the scripts folder
## of the package.
data(SpikeIn)
probes <- pm(SpikeIn) 
mas <- generateExprVal.method.mas(probes)
rl <- 2^robloxbioc(log2(t(probes)))
concentrations <- as.numeric(colnames(SpikeIn))
plot(concentrations, mas$exprs, log="xy", ylim=c(50,10000), type="b",
     ylab = "expression measures")
points(concentrations, rl[,1], pch = 20, col="orange", type="b")
legend("topleft", c("MAS", "roblox"), pch = c(1, 20))

\dontrun{
## "Not run" just because of computation time
require(affydata)
data(Dilution)
eset <- robloxbioc(Dilution)
## Affymetrix scale normalization
eset1 <- robloxbioc(Dilution, normalize = TRUE)
}

## using Illumina-Data
\dontrun{
## "Not run" just because of computation time
data(BLData)
BSData <- robloxbioc(BLData, eps.upper = 0.5)
}
}
\concept{normal location and scale}
\concept{infinitesimal robustness}
\concept{radius-minimax estimator}
\keyword{robust}
