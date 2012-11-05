##########################################
##                                      ## 
##        Tests for Expectation.R       ##
##                                      ##
##########################################

# .setUp(), .tearDown(): Either one or both functions have to be provided by the test case
#author, take precedence over the dummy definitions provided by the
#RUnit package and are called once for every test case identified.

 .setUp{ 

 
   ##expectation of Pareto distributed random variable
   expectation.Pareto = function(shape0=1,Min0=1){
    X = Pareto(shape=shape0,Min=Min0)
    return(E(X))  
   }

   test.expectationPareto = function(){
    checkEquals(expectation.Pareto(1,1), Inf)
    checkEquals(expectation.Pareto(2,1), 0)
   }    
  
#   test.HTMLInfo.Pareto = function(){
#    track <- tracker()
#    ## initialize the tracker
#    track$init()
#  
#    ## inspect the function
#    resFoo <- inspect(expectation.Pareto(1,1), track = track)
#    ## get the tracked function call info for all inspect calls
#    resTrack <- track$getTrackInfo()
#    }
# 
#  }

#  .tearDown(){
#   ##create HTML sites in folder ./results for all inspect calls
#   printHTML.trackInfo(resTrack,"TestSuite/TestExpectation")
#   }

}

#  
# Beispiele
# .setUp()
# {
# test.checkFunctions1 = function(){
# checkTrue(1 < 2, "check1") ## passes fine
# ## checkTrue(1 > 2, "check2") ## appears as failure in the test protocol
# v <- 1:3
# w <- 1:3
# checkEquals(v, w) ## passes fine
# names(v) <- c("A", "B", "C")
# ## checkEquals(v, w) ## fails because v and w have different names
# checkEqualsNumeric(v, w) ## passes fine because names are ignored
# x <- rep(1:12, 2)
# y <- rep(0:1, 12)
# res <- list(a=1:3, b=letters, LM=lm(y ~ x))
# res2 <- list(a=seq(1,3,by=1), b=letters, LM=lm(y ~ x))
# checkEquals(res, res2) ## passes fine
# checkIdentical(res, res)
# checkIdentical(res2, res2)
# ## checkIdentical(res, res2) ## fails because element ’a’ differs in type
# }
# }
# .tearDown()
# {}
# 
# fun <- function(x) {
# if(x)
# {
# stop("stop conditions signaled")
# }
# return()
# 
# }
# 
# .setUp()
# {
# test.checkFunctions2 = function(){
# checkException(fun(TRUE)) ## passes fine
# ## checkException(fun(FALSE)) ## failure, because fun raises no error
# checkException(fun(TRUE), silent=TRUE)
# ## special constants
# ## same behaviour as for underlying base functions
# checkEquals(NA, NA)
# checkEquals(NaN, NaN)
# checkEquals(Inf, Inf)
# checkIdentical(NA, NA)
# checkIdentical(NaN, NaN)
# checkIdentical(-Inf, -Inf)
# ## DEACTIVATED("here one can document on the reason for deactivation")
# }
# }
# 
# .tearDown()
# {}