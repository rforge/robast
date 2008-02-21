###############################################################################
## Some simple tests
###############################################################################

library(RobLox)

## sample
x <- rnorm(10000, mean = -2, sd = 3)

## location and scale, radius unknown
res1 <- roblox(x, returnIC = TRUE)
checkIC(res1$optIC)

res11 <- roblox(x, returnIC = TRUE, k = 2)
checkIC(res11$optIC)

res12 <- roblox(x, returnIC = TRUE, k = 5)
checkIC(res12$optIC)

roblox(x)
roblox(x, k = 2)
roblox(x, k = 5)


## location and scale, radius interval
res2 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(res2$optIC)

res21 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 2)
checkIC(res21$optIC)

res22 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 4)
checkIC(res22$optIC)

roblox(x, eps.lower = 0.15, eps.upper = 0.3)
roblox(x, eps.lower = 0.15, eps.upper = 0.3, k = 2)
roblox(x, eps.lower = 0.15, eps.upper = 0.3, k = 4)


## scale, radius interval
res3 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(res3$optIC)

res31 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 3)
checkIC(res31$optIC)

res32 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 6)
checkIC(res32$optIC)

roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3)
roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, k = 3)
roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, k = 6)


## location, radius interval
res4 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(res4$optIC)

res41 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 2)
checkIC(res41$optIC)

res42 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 5)
checkIC(res42$optIC)

roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3)
roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, k = 2)
roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, k = 5)


## some timings
system.time(for(i in 1:100) roblox(x, eps = 0.02))
system.time(for(i in 1:100) roblox(x))
system.time(for(i in 1:100) roblox(x, k = 2))
system.time(for(i in 1:100) roblox(x, k = 5))
system.time(for(i in 1:100) roblox(x, returnIC = TRUE))
