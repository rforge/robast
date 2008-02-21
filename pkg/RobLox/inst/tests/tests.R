###############################################################################
## Some simple tests
###############################################################################

library(RobLox)

x <- rnorm(10000, mean = -2, sd = 3)
res1 <- roblox(x, returnIC = TRUE)
checkIC(res1$optIC)

res11 <- roblox(x, returnIC = TRUE, k = 2)
checkIC(res11$optIC)

res12 <- roblox(x, returnIC = TRUE, k = 5)
checkIC(res12$optIC)
res1
res11
res12

res2 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(res2$optIC)

res21 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 2)
checkIC(res21$optIC)

res22 <- roblox(x, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 4)
checkIC(res22$optIC)
res2
res21
res22

res3 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(res3$optIC)

res31 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 3)
checkIC(res31$optIC)

res32 <- roblox(x, mean = -2, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 6)
checkIC(res32$optIC)
res3
res31
res32

res4 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE)
checkIC(res4$optIC)

res41 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 2)
checkIC(res41$optIC)

res42 <- roblox(x, sd = 3, eps.lower = 0.15, eps.upper = 0.3, returnIC = TRUE, k = 5)
checkIC(res42$optIC)
res4
res41
res42

system.time(for(i in 1:100) roblox(x, eps = 0.02))
system.time(for(i in 1:100) roblox(x))
system.time(for(i in 1:100) roblox(x, k = 2))
system.time(for(i in 1:100) roblox(x, k = 5))
system.time(for(i in 1:100) roblox(x, returnIC = TRUE))
