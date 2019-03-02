require(RobExtremes)
#
### plotIC
#
# GPD
fam <- GParetoFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y <- r(Y)(1000)
PlotIC(IC, y, withCall = FALSE)
#
# GEV
fam <- GEVFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y <- r(Y)(1000)
PlotIC(IC, y, withCall = FALSE)
#
# Weibull
fam <- WeibullFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y <- r(Y)(1000)
PlotIC(IC, y, withCall = FALSE)
#
### InfoPlot
#
# GPD
fam <- GParetoFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y <- r(Y)(1000)
InfoPlot(IC, y, withCall = FALSE)
#
# GEV
fam <- GEVFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y <- r(Y)(1000)
InfoPlot(IC, y, withCall = FALSE)
#
# Weibull
fam <- WeibullFamily()
IC <- optIC(model = fam, risk = asCov())
Y=distribution(fam)
y <- r(Y)(1000)
InfoPlot(IC, y, withCall = FALSE)
#
### ComparePlot
#
# GPD
fam <- GParetoFamily()
IC1 <- optIC(model = fam, risk = asCov())
IC2 <- makeIC(list(function(x)sin(x),function(x)x^2), L2Fam = fam)
Y=distribution(fam)
y <- r(Y)(1000)
ComparePlot(IC1, IC2, y, withCall = FALSE)
#
# GEV
fam <- GEVFamily()
IC1 <- optIC(model = fam, risk = asCov())
IC2 <- makeIC(list(function(x)sin(x),function(x)x^2), L2Fam = fam)
Y=distribution(fam)
y <- r(Y)(1000)
ComparePlot(IC1, IC2, y, withCall = FALSE)
#
# Weibull
fam <- WeibullFamily()
IC1 <- optIC(model = fam, risk = asCov())
IC2 <- makeIC(list(function(x)sin(x),function(x)x^2), L2Fam = fam)
Y=distribution(fam)
y <- r(Y)(1000)
ComparePlot(IC1, IC2, y, withCall = FALSE)
#
### CniperPointPlot
#
# GPD
L2fam <- GParetoFamily()
CniperPointPlot(fam=L2fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)
#
# GEV
L2fam <- GEVFamily()
CniperPointPlot(fam=L2fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)
#
# Weibull
L2fam <- WeibullFamily()
CniperPointPlot(fam=L2fam, main = "Gamma", lower = 0, upper = 5, withCall = FALSE)
