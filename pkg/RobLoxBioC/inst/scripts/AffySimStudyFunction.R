###############################################################################
## Function to perform simulation study comparing Tukey's biweight with the 
## rmx estimator
###############################################################################

## n: sample size
## M: Monte Carlo replications
## eps: amount of contamination
## seed: seed for simulations
## eps.lower: eps.lower for rmx estimator
## eps.upper: eps.upper for rmx estimator
## steps: number of steps used for estimator construction
## fsCor: perform finite-sample correction
## contD: contaminating distribution
## plot1: plot densities of ideal and real situation
## plot2: plot 20 randomly chosen samples out of the M samples
## plot3: boxplots of the estimates
AffySimStudy <- function(n, M, eps, seed = 123, eps.lower = 0, eps.upper = 0.2, 
                         steps = 3, fsCor = TRUE, contD, 
                         plot1 = TRUE, plot2 = FALSE, plot3 = TRUE){
    if(plot1){
        from <- min(-6, q(contD)(1e-15))
        to <- max(6, q(contD)(1-1e-15))
        curve(dnorm, from = from, to = to, lwd = 2, n = 201, 
              main = "Comparison: ideal vs. real")
        fun <- function(x) (1-eps)*dnorm(x)+eps*d(contD)(x)
        curve(fun, from = from, to = to, add = TRUE, col = "orange", 
              lwd = 2, n = 201)
        legend("topleft", legend = c("ideal", "real"), 
              fill = c("black", "orange"))
    }

    set.seed(seed)
    r <- rbinom(n*M, prob = eps, size = 1)
    Mid <- rnorm(n*M)
    Mcont <- r(contD)(n*M)
    Mre <- matrix((1-r)*Mid + r*Mcont, ncol = n)

    if(plot2){
        library(lattice)
        ind <- sample(1:M, 20)
        if(plot1) dev.new()
        print(
          stripplot(rep(1:20, each = 20) ~ as.vector(Mre[ind,]), 
          ylab = "samples", xlab = "x", pch = 20,
          main = "Randomly chosen samples")
        )
    }

    ## ML-estimator: mean and sd
    Mean <- rowMeans(Mre)
    Sd <- sqrt(rowMeans((Mre-Mean)^2))
    ## Median and MAD
    Median <- rowMedians(Mre)
    Mad <- rowMedians(abs(Mre - Median))/qnorm(0.75)
    ## Tukey 1-step + MAD
    Tukey <- apply(Mre, 1, function(x) tukey.biweight(x))
    Tukey <- cbind(Tukey, Mad)

    ## Radius-minimax estimator
    RadMinmax1 <- estimate(rowRoblox(Mre, eps.lower = eps.lower, 
                                    eps.upper = eps.upper, k = steps,
                                    fsCor = fsCor))
    RadMinmax2 <- estimate(rowRoblox(Mre, sd = Mad, eps.lower = eps.lower, 
                                    eps.upper = eps.upper, k = steps,
                                    fsCor = fsCor))
    RadMinmax2 <- cbind(RadMinmax2, Mad)

    if(plot3){
        Ergebnis1 <- list(Mean, Median, Tukey[,1], RadMinmax1[,1], RadMinmax2[,1])
        Ergebnis2 <- list(Sd, Mad, Tukey[,2], RadMinmax1[,2], RadMinmax2[,2])
        myCol <- brewer.pal(4, "Dark2")
        if(plot1 || plot2) dev.new()
        layout(matrix(c(1, 1, 1, 1, 3, 2, 2, 2, 2, 3), ncol = 2))
        boxplot(Ergebnis1, col = myCol, pch = 20, main = "Location")
        abline(h = 0)
        boxplot(Ergebnis2, col = myCol, pch = 20, main = "Scale")
        abline(h = 1)
        op <- par(mar = rep(2, 4))
        plot(c(0,1), c(1, 0), type = "n", axes = FALSE)
        legend("center", c("ML", "Med/MAD", "biweight", "rmx", "rmx/MAD"),
               fill = myCol, ncol = 5, cex = 1.5)
        par(op)
    }

    ## ML-estimator
    MSE1.1 <- n*mean(Mean^2)
    ## Median + MAD
    MSE2.1 <- n*mean(Median^2)
    ## Tukey
    MSE3.1 <- n*mean(Tukey[,1]^2)
    ## Radius-minimax
    MSE4.1 <- n*mean(RadMinmax1[,1]^2)
    MSE5.1 <- n*mean(RadMinmax2[,1]^2)
    empMSE <- data.frame(ML = MSE1.1, Med = MSE2.1, Tukey = MSE3.1, 
                         "rmx" = MSE4.1, "rmx1" = MSE5.1)
    rownames(empMSE) <- "n x empMSE (loc)"

    ## ML-estimator
    MSE1.2 <- n*mean((Sd-1)^2)
    ## Median + MAD
    MSE2.2 <- n*mean((Mad-1)^2)
    ## Tukey
    MSE3.2 <- MSE2.2
    ## Radius-minimax
    MSE4.2 <- n*mean((RadMinmax1[,2]-1)^2)
    MSE5.2 <- n*mean((RadMinmax2[,2]-1)^2)
    empMSE <- rbind(empMSE, c(MSE1.2, MSE2.2, MSE3.2, MSE4.2, MSE5.2))
    rownames(empMSE)[2] <- "n x empMSE (scale)"
    empMSE <- rbind(empMSE, c(MSE1.1 + MSE1.2, MSE2.1 + MSE2.2, MSE3.1 + MSE3.2, 
                              MSE4.1 + MSE4.2, MSE5.1 + MSE5.2))
    rownames(empMSE)[3] <- "n x empMSE (loc + scale)"

    empMSE
}


