.getTau <- function(data, model, level, rob=TRUE, of.interest, what){
  if(rob){
     est <- roptest(data,model, risk=RMXRRisk())
     L2FamC <- CallL2Fam(pIC(est))
  }else{
     est <- MLEstimator(data,model)
     IC <- optIC(model,risk=asCov())
     L2FamC <- CallL2Fam(IC)
     L2FamC$scale <- estimate(est)["scale"]
     L2FamC$shape <- estimate(est)["shape"]
  }
  eval(what)
  L2FamC$of.interest <- of.interest # "quantile"
  L2Fam <- eval(L2FamC)
  res <- param(L2Fam)@trafo(estimate(est))
  VaR <- res[[1]]
  varVaR <- (res[[2]]) %*% asvar(est) %*% t(res[[2]])
  res <- c(VaR,sqrt(varVaR/length(data)))
  names(res) <- c("Risk","varofRisk")
  class(res) <- "riskMeasure"
  res
}
print.riskMeasure <- function(x, level=NULL, ...){
   mc <- as.list(match.call(expand.dots=TRUE)[-1])
   digits <- if(is.null(mc$digits)) 3 else  mc$digits
   if(is.null(level)){
      cat(" ",signif(x[1],digits),"\n")
      cat("(",signif(x[2],digits),")\n")
   }else{qn <- qnorm((level+1)/2)
      CI <- c(-1,1)*qn*x[2]+x[1]
      cat(" ",signif(x[1],digits),"         [", signif(CI[1],digits), ",",
              signif(CI[2],digits),"]\n")
  }
}


getVaR <- function(data, model, level, rob=TRUE)
             .getTau(data, model, level, rob, of.interest="quantile", substitute(L2FamC$p <- level))

getCVaR <- function(data, model, level, rob=TRUE)
             .getTau(data, model, level, rob, of.interest="expected shortfall", substitute(L2FamC$p <- level))

getEL <- function(data, model, N0, rob=TRUE)
             .getTau(data, model, N0, rob, of.interest="expected loss", substitute(L2FamC$N <- N0))


