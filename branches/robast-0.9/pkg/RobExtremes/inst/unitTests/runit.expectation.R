##########################################
##                                      ## 
##        Tests for Expectation.R       ##
##                                      ##
##########################################

# .setUp(), .tearDown():
# Either one or both functions have to be provided by the test case
# author, take precedence over the dummy definitions provided by the
# RUnit package and are called once for every test case identified.

# we construct different objects for testing the expectation operator
.setUp <- function() {
  # expectation of Pareto distributed random variable
  expectation.Pareto <<- function(shape0=1,Min0=1){
    X = Pareto(shape=shape0,Min=Min0)
    return(E(X))  
  }
}

# test for the expectation of the pareto-distribution
test.expectationPareto1 <- function() {
  checkEquals(expectation.Pareto(1, 1), Inf)
}

test.expectationPareto2 <- function() {
  checkEquals(expectation.Pareto(2, 1), 0)
}
