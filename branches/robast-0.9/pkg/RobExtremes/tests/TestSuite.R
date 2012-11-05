##########################################
##                                      ## 
## TestSuite for RobExtremes Package    ##
##                                      ##
##########################################

##load RobExtremes package
require(RobExtremes)

##load RUnit package
require(RUnit)

#Options: if the exceptions are shown as errors
## Not run:
## quiet log output
ro <- getOption("RUnit")
ro$silent <- TRUE
ro$verbose <- 0L
options("RUnit"=ro)
## End(Not run)

## run test suite
myTestSuite <- defineTestSuite("RUnit Tests for RobExtremes",dirs="./TestSuite",
testFileRegexp = "^Test.+R$",testFuncRegexp = "^test.+")
testResult <- runTestSuite(myTestSuite)

## is test suite valid?
isValidTestSuite(myTestSuite)

## get all the errors
#getErrors(myTestSuite)

## prints detailed text protocol
## to standard out:
printTextProtocol(testResult, "testResult.txt",showDetails = TRUE)
printHTMLProtocol(testResult, "testResult.html")


#  
# runTestFile(absFileName, useOwnErrorHandler = TRUE,
#                  testFuncRegexp = "^test.+",
#                  rngKind = "Marsaglia-Multicarry",
#                  rngNormalKind = "Kinderman-Ramage",
#                  verbose = getOption("RUnit")$verbose)

