#' # GPD
#' fam = GParetoFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' data = r(Y)(1000)
#' InfoPlot(IC, data, withCall = FALSE)
#'
#' # GEV
#' fam = GEVFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' data = r(Y)(1000)
#' InfoPlot(IC, data, rescale = TRUE, withCall = FALSE)
#'
#' # Weibull
#' fam = WeibullFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' data = r(Y)(1000)
#' InfoPlot(IC, data, withCall = FALSE)
#' # GPD
#' fam = GParetoFamily()
#' CniperPointPlot(fam=fam, main = "GPD", lower = 0, upper = 10, withCall = FALSE)
#' # GEV
#' fam = GEVFamily()
#' CniperPointPlot(fam=fam, main = "GEV", lower = 0, upper = 5, withCall = FALSE)
#' # Gamma
#' fam = GammaFamily()
#' CniperPointPlot(fam=fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)
#' # Weibull
#' fam = WeibullFamily()
#' CniperPointPlot(fam=fam, main = "Weibull", withCall = FALSE)
#' @examples
#' # GPD
#' fam = GParetoFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' PlotIC(IC, y, withCall = FALSE)
#'
#' # GEV
#' fam = GEVFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' PlotIC(IC, y, rescale = TRUE, withCall = FALSE)
#'
#' # Gamma
#' fam = GammaFamily()
#' rfam = InfRobModel(fam, ContNeighborhood(0.5))
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' PlotIC(IC, y, withCall = FALSE)
#'
#' # Weibull
#' fam = WeibullFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' PlotIC(IC, y, withCall = FALSE)
#' @examples
#' # GPD
#' fam = GParetoFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' ComparePlot(IC, y, withCall = FALSE)
#'
#' # GEV
#' fam = GEVFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' ComparePlot(IC, y, rescale = TRUE, withCall = FALSE)
#'
#' # Gamma
#' fam = GammaFamily()
#' rfam = InfRobModel(fam, ContNeighborhood(0.5))
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' ComparePlot(IC, y, withCall = FALSE)
#'
#' # Weibull
#' fam = WeibullFamily()
#' IC <- optIC(model = fam, risk = asCov())
#' Y=distribution(fam)
#' y = r(Y)(1000)
#' ComparePlot(IC, y, withCall = FALSE)
